# Leopard2 current-literature refresh

## Scope and result

Search date: 2026-07-16 UTC.

This is the R23 refresh required before design freeze.  It searched for work on
Reed-Solomon erasure coding, LCH-FFT/additive transforms, truncated transforms,
and later papers related to R10-R22.  It is a bounded database and primary-source
review, not a systematic-review protocol and not evidence that Leopard2 or any
candidate is state of the art.

No reviewed paper requires changing the production high/low profile definitions,
the legacy compatibility gate, or the specialized R10 decoder plan.  One material
source omission was found: the 2023 ISIT paper on Cantor-constructed fields and
truncated LCH transforms should be read in full before promoting experiments C1-C3
or claiming an O(1)-workspace/message-only implementation.  The 2024 rateless and
four-parity papers add useful experimental comparisons, but their coordinate,
field-extension, or generator constructions are not silently compatible with a
Leopard wire profile.

## Exact endpoint queries and access results

The following URLs and queries were issued on 2026-07-16 UTC.

### Required R23 endpoints

1. arXiv

   - URL: https://arxiv.org/search/?query=Reed-Solomon+LCH-FFT&searchtype=all
   - Exact query: `Reed-Solomon LCH-FFT`, all fields.
   - Result: the endpoint was readable and reported no results for that literal
     tokenized query.  This does not show that arXiv contains no related papers;
     several relevant papers use different terminology or are publisher/author
     copies rather than arXiv submissions.
   - Follow-up URLs:
     - https://arxiv.org/search/?query=%22Reed-Solomon%22+%22erasure%22+%22FFT%22&searchtype=all&abstracts=show&order=-announced_date_first&size=50
     - https://arxiv.org/search/?query=%22Reed-Solomon%22+%22additive+FFT%22&searchtype=all&abstracts=show&order=-announced_date_first&size=50
     - https://arxiv.org/search/?query=%22Reed-Solomon%22+%22truncated%22+transform&searchtype=all&abstracts=show&order=-announced_date_first&size=50
   - Result: all three follow-up search fetches timed out in the current web
     environment.  Attempts to use `export.arxiv.org/api/query` were rejected by
     the browsing safety layer.  Known arXiv records R13, R15-R17, and R19 remain
     in the main bibliography; the failed refresh endpoints are not evidence
     against them.

2. IEEE Xplore

   - URL: https://ieeexplore.ieee.org/search/searchresult.jsp?newsearch=true&queryText=Reed-Solomon%20LCH-FFT
   - Exact query: `Reed-Solomon LCH-FFT`.
   - Result: the aggregate search endpoint returned an internal fetch error in
     this environment.  Individual official IEEE records found through exact-title
     queries were reachable, including documents 10496475 and 10827202.  Some
     records exposed metadata only, without accessible full text.
   - Follow-up exact-title queries:
     - `"Two Fast Erasure Decoding Algorithms for Reed-Solomon Codes Based on LCH-FFT"`
     - `"Reduced-Complexity Erasure Decoding" "LCH-FFT"`
     - `"An Efficient Reed-Solomon Erasure Code over Cantor-constructed Binary Extension Finite Fields"`
     - `"Fast Erasure Decoding of Reed-Solomon Codes Based on Decomposition of Vandermonde Matrix"`
     - `"Construction of Reed-Solomon Erasure Codes with Four Parities"`

3. DBLP

   - URL: https://dblp.org/search?q=Reed-Solomon%20LCH-FFT
   - Exact query: `Reed-Solomon LCH-FFT`.
   - Result: direct aggregate fetch returned an internal error in this
     environment.  Search-index and individual DBLP records were reachable.  The
     exact term located the final R10 journal record, its 2023 low-rate precursor,
     and the foundational R16-R17 records.  Attempts to open DBLP JSON search API
     URLs were rejected by the browsing safety layer.
   - Individual records used:
     - https://dblp.org/rec/journals/tit/ChenLTHCYLBB26
     - https://dblp.org/rec/conf/isit/ChenLLCHB23
     - https://dblp.org/rec/conf/isit/LiHLC23
     - https://dblp.org/rec/journals/tcom/YuLH24
     - https://dblp.org/rec/conf/wcsp/LinLLSLH24
     - https://dblp.org/rec/conf/isit/ChenLHB25

### Supplemental exact web queries

The search was broadened with these exact query strings because the aggregate
endpoints were incomplete or inaccessible:

- `site:arxiv.org Reed-Solomon "LCH-FFT" 2024 2025 2026`
- `site:ieeexplore.ieee.org/document Reed-Solomon "LCH-FFT"`
- `site:dblp.org "LCH-FFT" "Reed-Solomon"`
- `"Two Fast Erasure Decoding Algorithms for Reed-Solomon Codes Based on LCH-FFT" citations`
- `"Reed-Solomon" "LCH-FFT" 2025`
- `"Reed-Solomon" "LCH-FFT" 2026`
- `"improved truncated LCH transform"`
- `"A Class of Rateless Reed-Solomon Codes with Near-Linear Computational Complexities"`
- `"Construction of Reed-Solomon Erasure Codes with Four Parities"`
- `"Efficient Decoding of a Class of Reed-Solomon Codes over Fermat Fields"`
- `"A Fast Chinese Remaindering Transform Over Finite Fields"`

Candidate metadata was then checked against DOI landing pages, DBLP records,
official journal pages, or author-hosted manuscripts where available.  Search
snippets and ResearchGate were used only to discover titles or inaccessible
records, not as the final authority for production mathematics.

## Candidates reviewed

### Final R10 publication

Title: *Two Fast Erasure Decoding Algorithms for Reed-Solomon Codes Based on
LCH-FFT*.

Authors: Chao Chen, Sian-Jheng Lin, Nianqi Tang, Yunghsiang S. Han, Suihua Cai,
Leilei Yu, Zhongwei Li, Baoming Bai, and Bo Bai.

Year and venue: 2026, IEEE Transactions on Information Theory 72(6), 3784-3798.

URLs:

- DOI: https://doi.org/10.1109/TIT.2026.3685291
- DBLP: https://dblp.org/rec/journals/tit/ChenLTHCYLBB26
- Author manuscript: https://i4ai.org/hanyunghsiang/IT2026.pdf

Relevance: this is R10 itself, now with final journal metadata.  It remains the
primary source for the O(N log K) low-rate and O(N log(N-K)) high-rate erasure
decoders and their message-only specialization.

Plan decision: no algorithm change.  Cite the final DOI/volume/pages in release
documentation while retaining the author manuscript URL used for the full-text
derivation.  Because the article appeared in June 2026, absence of a later
indexed improvement on 2026-07-16 has little evidentiary weight.

### Reduced-complexity low-rate precursor

Title: *Reduced-Complexity Erasure Decoding of Low-Rate Reed-Solomon Codes Based
on LCH-FFT*.

Authors: Chao Chen, Sian-Jheng Lin, Zhongwei Li, Suihua Cai, Yunghsiang S. Han,
and Bo Bai.

Year and venue: 2023, IEEE ISIT, 1015-1019.

URLs:

- DOI: https://doi.org/10.1109/ISIT54713.2023.10206549
- DBLP: https://dblp.org/rec/conf/isit/ChenLLCHB23

Relevance: conference precursor for the low-rate branch later developed in R10.

Plan decision: no change; R10 is the fuller and newer production source.  Keep
the precursor in the provenance trail rather than treating it as an independent
algorithm candidate.

### Cantor-field truncated LCH erasure code

Title: *An Efficient Reed-Solomon Erasure Code over Cantor-constructed Binary
Extension Finite Fields*.

Authors: Zhongwei Li, Yunghsiang Sam Han, Sian-Jheng Lin, and Chao Chen.

Year and venue: 2023, IEEE ISIT, 826-831.

URLs:

- DOI: https://doi.org/10.1109/ISIT54713.2023.10206562
- DBLP: https://dblp.org/rec/conf/isit/LiHLC23

Relevance: bibliographic and abstract records describe an improved truncated LCH
transform for discrete intervals, O(n log T) erasure-code work, constant
auxiliary-space claims, and a message-only variation.  Those topics directly
overlap C1-C3, dependency pruning, tiled scratch, and original-only recovery.
The full paper was not accessible from the reviewed official endpoints, so no
formula from it has been accepted on abstract evidence.

Plan decision: no frozen production-profile change, but this changes source
priority.  Obtain and read the primary full text before promoting C1-C3 or making
space-complexity claims.  Any coordinate construction must still pass the
byte-for-byte legacy wire test; a Cantor-field relationship alone is not proof
of generator equivalence.

### Rateless LCH Reed-Solomon codes

Title: *A Class of Rateless Reed-Solomon Codes With Near-Linear Computational
Complexities*.

Authors: Leilei Yu, Sian-Jheng Lin, and Yunghsiang S. Han.

Year and venue: 2024, IEEE Transactions on Communications 72(11), 6677-6687.

URLs:

- DOI: https://doi.org/10.1109/TCOMM.2024.3405336
- DBLP: https://dblp.org/rec/journals/tcom/YuLH24
- Author manuscript: https://i4ai.org/hanyunghsiang/tcom2024-3.pdf

Relevance: the primary manuscript describes Cantor field-tower extensions and an
LCH schedule that produces encoded packets on demand instead of only in
power-of-two bursts.  This is relevant to progressive encoding and field-boundary
research.

Plan decision: add it to the source set for Experiment I and future rateless work,
but do not alter legacy high or low V1.  Preservative field extension and
on-demand coordinate generation define a different persistent code identity and
would require a separately versioned profile, interoperability specification,
and MDS tests.

### Four-parity systematic Vandermonde construction

Title: *Construction of Reed-Solomon Erasure Codes With Four Parities Based on
Systematic Vandermonde Matrices*.

Authors: Leilei Yu and Yunghsiang S. Han.

Year and venue: 2024, IEEE Transactions on Computers 73(7), 1875-1882.

URLs:

- DOI: https://doi.org/10.1109/TC.2024.3387069
- IEEE: https://ieeexplore.ieee.org/document/10496475/
- Author manuscript: https://i4ai.org/hanyunghsiang/TC2024.pdf

Relevance: an algebraic, four-parity systematic Vandermonde/RM-transform
construction.  It is pertinent to tiny-R direct kernels and is a focused follow-up
to the RM-transform family represented by R18.

Plan decision: include as a comparison for R=4 and Experiment O.  Do not import
its generator under legacy high V1 without a generator-matrix equivalence proof;
otherwise it is a new profile.  It does not replace the general direct LxL repair
path or the R10 high-rate decoder.

### Vandermonde-decomposition erasure decoder

Title: *Fast Erasure Decoding of Reed-Solomon Codes Based on Decomposition of
Vandermonde Matrix*.

Authors: Zebing Lin, Jingjie Lv, Pingping Li, Linqi Song, Hui Liang, and Hanxu Hou.

Year and venue: 2024, 16th International Conference on Wireless Communications
and Signal Processing, 96-101.

URLs:

- DOI: https://doi.org/10.1109/WCSP62071.2024.10827202
- IEEE: https://ieeexplore.ieee.org/document/10827202
- DBLP: https://dblp.org/rec/conf/wcsp/LinLLSLH24

Relevance: the accessible metadata identifies a Vandermonde-decomposition/LU
erasure decoder and places it in the RM-transform comparison family.  The full
text was not accessible through the reviewed official endpoint, so detailed
complexity or wire claims were not adopted.

Plan decision: optional comparison cell for direct/small-loss and R18 experiments.
No production change without the paper, an independent derivation, and a fair
same-code benchmark.

### LCH-FFT algorithms and FPGA architectures

English title: *New Algorithms and Architectures for Reed-Solomon Codes Based on
LCH-FFT*.

Authors: Chao Chen, Nianqi Tang, Yunghsiang Han, Xiaotian Wang, and Baoming Bai.

Year and venue: 2025, Journal of Tsinghua University (Science and Technology)
65(11), 2053-2066.

URL: https://jst.tsinghuajournals.com/CN/10.16511/j.cnki.qhdxxb.2025.27.045

Relevance: the official bilingual abstract describes two systematic encoder
architectures, an FDMA-based decoder architecture, conversion for cyclic RS
codes, and an FPGA RS(544,514) evaluation.  The decoder includes unknown-error
machinery and the evidence is hardware-oriented.

Plan decision: useful for encoder architecture and later hardware experiments.
It does not move FDMA or key-equation work into the production erasure hot path,
and it does not establish a software or legacy-wire improvement.

### Fast Chinese remaindering transform

Title: *A Fast Chinese Remaindering Transform Over Finite Fields*.

Authors: Chao Chen, Sian-Jheng Lin, Yunghsiang S. Han, and Baoming Bai.

Year and venue: 2025, IEEE ISIT, 1-6.

URLs:

- DOI: https://doi.org/10.1109/ISIT63088.2025.11195593
- DBLP: https://dblp.org/rec/conf/isit/ChenLHB25

Relevance: a newer finite-field transform by authors in the LCH lineage.  Only
bibliographic metadata was accessible; no primary text reviewed here connected
it to the exact Leopard erasure-code maps.

Plan decision: monitor for the exact-size transform experiments, but no plan
change on title/metadata alone.

### Reviewed but out of the production erasure scope

*Efficient Decoding of a Class of Reed-Solomon Codes over Fermat Fields* by Chao
Chen, Baoming Bai, Xiao Ma, Yunghsiang S. Han, Nianqi Tang, and Xiaotian Wang,
IEEE ISIT 2024, DOI https://doi.org/10.1109/ISIT57864.2024.10619665, concerns a
different field family, syndrome/Chien-search acceleration, and error-decoder
architecture.  It does not change the GF8/GF16 LCH erasure plan.

R11 (Parallel Welch-Berlekamp, 2025), R12 (fast error-and-erasure decoding,
2024), and R19 (generalized RS/alternant decoding, 2025) were also rediscovered
through author/DBLP searches.  They are already in the main bibliography and
remain outside the normal pure-erasure hot path as originally scoped.

## Design-freeze decision

The 2026-07-16 refresh makes these bounded adjustments:

1. Use final R10 journal metadata in citations; keep R10 as the specialized
   decoder authority.
2. Add Li-Han-Lin-Chen ISIT 2023 to the required reading gate for C1-C3 and
   low-memory/message-only implementation claims.
3. Add Yu-Lin-Han 2024 to Experiment I/future separately versioned rateless work.
4. Add Yu-Han 2024 and Lin-Lv-Li-Song-Liang-Hou 2024 as R=4/direct/RM-transform
   comparison candidates.
5. Keep the 2025 architecture and Chinese-remaindering papers experimental until
   primary mathematical text, software-equivalent profiles, and target-machine
   measurements justify more.

Nothing found authorizes a wire change, removal of the direct scalar oracle,
replacement of shortening/puncturing, or integration of unknown-error decoding
into erasure execution.  Because several aggregate endpoints failed and R10 was
only recently published, repeat this refresh before a public release.
