
import time
import requests
import json

def _extract_sublocs_from_comments(comments):
    """Try many shapes of UniProt 'comments' entries to extract subcellular locations."""
    found = []
    for c in comments or []:
        # commentType might be "SUBCELLULAR LOCATION" or variants
        ctype = c.get("commentType") or c.get("type") or ""
        if isinstance(ctype, str) and "subcellular" in ctype.lower():
            # 1) modern structured form: 'subcellularLocations' (list)
            sitems = c.get("subcellularLocations") or c.get("subcellularLocation") or c.get("subcellular_locations") or []
            if isinstance(sitems, dict):
                sitems = [sitems]
            for s in sitems or []:
                # location may be nested under 'location' or 'locationName' or 'location' dict with 'value'
                loc = None
                if isinstance(s, dict):
                    # common path: s['location']['value']
                    loc_obj = s.get("location") or s.get("locationName") or s.get("locationText") or {}
                    if isinstance(loc_obj, dict):
                        loc = loc_obj.get("value") or loc_obj.get("name") or loc_obj.get("text")
                    elif isinstance(loc_obj, str):
                        loc = loc_obj
                    # some entries include qualifiers like 'note' or 'topology' — we can include them if needed
                    # also location may be direct string fields
                    if not loc:
                        # try keys directly
                        for k in ("value", "name", "text", "location"):
                            if isinstance(s.get(k), str) and s.get(k).strip():
                                loc = s.get(k).strip()
                                break
                elif isinstance(s, str):
                    loc = s
                if loc:
                    found.append(str(loc).strip())
            # 2) fallback: 'texts' or 'text' fields containing human-readable description
            if not found:
                texts = c.get("texts") or []
                if isinstance(texts, dict):
                    texts = [texts]
                for t in texts or []:
                    if isinstance(t, dict):
                        tv = t.get("value") or t.get("text") or t.get("location")
                        if tv:
                            found.append(str(tv).strip())
                    elif isinstance(t, str):
                        found.append(t.strip())
            # 3) sometimes the comment has 'note' or 'subcellularLocation' directly
            if not found:
                alt = c.get("subcellularLocation") or c.get("note")
                if isinstance(alt, str) and alt.strip():
                    found.append(alt.strip())
    # dedupe preserving order
    seen = []
    for s in found:
        s2 = s.strip()
        if s2 and s2 not in seen:
            seen.append(s2)
    return seen

def _extract_sublocs_from_features(features):
    """Some UniProt entries may use features to denote subcellular location."""
    found = []
    for f in features or []:
        # features have 'type' like 'subcellular location' in some variants
        ftype = f.get("type") or ""
        if isinstance(ftype, str) and "subcellular" in ftype.lower():
            # look for 'description', 'location', 'note'
            for k in ("description", "location", "note", "text"):
                val = f.get(k)
                if isinstance(val, str) and val.strip():
                    found.append(val.strip())
                    break
            # sometimes location is nested
            loc = f.get("location")
            if isinstance(loc, dict):
                for k in ("value","name","text"):
                    if loc.get(k):
                        found.append(str(loc.get(k)).strip())
                        break
    # dedupe
    seen = []
    for s in found:
        if s not in seen:
            seen.append(s)
    return seen

def extract_subcellular_from_json_more(data):
    """
    Try many JSON paths to extract subcellular location strings.
    Returns semicolon-separated string or None.
    """
    # 1) comments variant
    comments = data.get("comments") or data.get("comment") or []
    locs = _extract_sublocs_from_comments(comments)
    # 2) features variant
    if not locs:
        features = data.get("features") or []
        locs = _extract_sublocs_from_features(features)
    # 3) top-level 'subcellularLocation' field (less common but possible)
    if not locs:
        top = data.get("subcellularLocation") or data.get("subcellular_location")
        if top:
            if isinstance(top, list):
                for t in top:
                    if isinstance(t, dict):
                        loc = t.get("location") or t.get("locationName") or t.get("value") or t.get("text")
                        if isinstance(loc, dict):
                            for k in ("value","name","text"):
                                if loc.get(k):
                                    locs.append(str(loc.get(k)).strip())
                                    break
                        elif isinstance(loc, str) and loc.strip():
                            locs.append(loc.strip())
                    elif isinstance(t, str) and t.strip():
                        locs.append(t.strip())
            elif isinstance(top, dict):
                loc = top.get("location") or top.get("value") or top.get("text")
                if isinstance(loc, dict):
                    for k in ("value","name","text"):
                        if loc.get(k):
                            locs.append(str(loc.get(k)).strip())
                            break
                elif isinstance(loc, str) and loc.strip():
                    locs.append(loc.strip())
            elif isinstance(top, str) and top.strip():
                locs.append(top.strip())
    # 4) fallback: check 'comments' text blobs more widely (search for key phrases)
    if not locs:
        for c in comments or []:
            # some older formats put location in c['text'][0]['value'] or similar
            texts = c.get("texts") or c.get("text") or []
            if isinstance(texts, dict):
                texts = [texts]
            for t in texts or []:
                if isinstance(t, dict):
                    val = t.get("value") or t.get("text")
                    if val and isinstance(val, str):
                        # try to extract short location phrases heuristically: look for known organelles
                        organelles = ["mitochondr", "nucle", "cytoplasm", "membran", "plasma membrane", "secreted", "extracellular", "endoplasmic", "golgi", "lysosom", "peroxisom"]
                        low = val.lower()
                        hits = []
                        for o in organelles:
                            if o in low:
                                # try to find the full word around it (rough)
                                hits.append(o)
                        if hits:
                            # add the whole text as fallback
                            locs.append(val.strip())
    # final tidy
    seen = []
    for s in locs:
        s2 = str(s).strip()
        if s2 and s2 not in seen:
            seen.append(s2)
    if not seen:
        return None
    # normalize short forms: attempt to split composite annotations
    # many entries contain "Mitochondrion; Mitochondrial inner membrane" etc already; preserve as-is
    return "; ".join(seen)

def lookup_uniprot_annotation_improved(accession, max_retries=4, timeout=30):
    """
    Per-accession JSON query; returns (gene_name_or_None, subcellular_location_or_None).
    More defensive and logs when no info present.
    """
    url = UNIPROT_ENDPOINT + f"{accession}.json"
    backoff = 1.0
    for attempt in range(1, max_retries+1):
        try:
            r = requests.get(url, timeout=timeout)
            # network down or blocked -> r may fail
            if r.status_code == 200:
                try:
                    data = r.json()
                except Exception:
                    return (None, None)
                # gene name
                gene_name = None
                genes = data.get("genes") or []
                if genes:
                    for g in genes:
                        if isinstance(g, dict):
                            gn = g.get("geneName") or g.get("gene_name") or {}
                            if isinstance(gn, dict):
                                val = gn.get("value") or gn.get("primary") or gn.get("name")
                                if val:
                                    gene_name = str(val)
                                    break
                            # synonyms
                            syns = g.get("synonyms") or []
                            if syns:
                                first = syns[0]
                                if isinstance(first, str):
                                    gene_name = first
                                    break
                                if isinstance(first, dict):
                                    val = first.get("value") or first.get("name")
                                    if val:
                                        gene_name = val
                                        break
                # subcellular location extraction (robust)
                subloc = extract_subcellular_from_json_more(data)
                return (gene_name, subloc)
            elif r.status_code == 404:
                return (None, None)
            elif r.status_code in (429, 503):
                time.sleep(backoff)
                backoff *= 1.7
            else:
                # other server error -> retry
                time.sleep(backoff)
                backoff *= 1.7
        except requests.RequestException:
            # network error -> retry after wait
            time.sleep(backoff)
            backoff *= 1.7
    return (None, None)
