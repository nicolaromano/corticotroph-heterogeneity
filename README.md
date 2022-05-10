This repository contains scripts for the analysis of corticotrophs in scRNAseq pituitary datasets.

The following datasets have been used

| PMID                                                 | Species | Reference                    | Data source                                                                           |
|------------------------------------------------------|---------|------------------------------|---------------------------------------------------------------------------------------|
| [33571131](https://pubmed.ncbi.nlm.nih.gov/33571131) |  Mouse  | Lopez et al.2021             |                                                                                       |
| [33373440](https://pubmed.ncbi.nlm.nih.gov/33373440) |  Mouse  | Allensworth-James et al.2021 |                                                                                       |
| [33808370](https://pubmed.ncbi.nlm.nih.gov/33808370) |  Mouse  | Scagliotti et al.2021        |                                                                                       |
| [34161279](https://pubmed.ncbi.nlm.nih.gov/34161279) |  Mouse  | Vennekens et al.2021         |                                                                                       |
| [31915267](https://pubmed.ncbi.nlm.nih.gov/31915267) |  Mouse  | Chen et al.2020              |                                                                                       |
| [32193873](https://pubmed.ncbi.nlm.nih.gov/32193873) |  Mouse  | Ho et al.2020                |                                                                                       |
| [31444346](https://pubmed.ncbi.nlm.nih.gov/31444346) |  Mouse  | Mayran et al.2019            | [PRJNA517137](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA517137&o=acc_s%3Aa) |
| [30335147](https://pubmed.ncbi.nlm.nih.gov/30335147) |  Mouse  | Cheung et al.2018            |                                                                                       |
| [35058881](https://pubmed.ncbi.nlm.nih.gov/35058881) |   Rat   | Kučka et al.2021             |                                                                                       |
| [31620083](https://pubmed.ncbi.nlm.nih.gov/31620083) |   Rat   | Fletcher et al.2019          | [PRJNA546549](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA546549&o=acc_s%3Aa) |

The scripts have been numbered sequentially to help navigating the repository. A brief description of each script follows

---

`00_download_data.sh`

Downloads raw data for each dataset - Note that this takes **a lot** of space.
Usage

    00_download_data.sh SRA_id output_dir
    
For example

    00_download_data.sh SRR9203724 /home/pitdatasets/Fletcher2019
