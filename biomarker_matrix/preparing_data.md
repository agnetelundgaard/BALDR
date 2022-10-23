# Fetching and preparing data from public sources
While many of these data sources offer API access, we instead opted for an approach where we downloaded bulk data and processed the bulk data locally. 
For some of the data sources, the process of obtaining and preparing this bulk data is very straightforward, whereas other sources require more handiwork.

All of these steps were performed on a Macbook Pro from 2017 with 16 GB of memory, using MacOS Monterey.
More detailed information on R and package versions can be found at the bottom of this report.

The instructions should be valid for the most recent data versions as of August 5th, 2022.

### Uniprot

##### Data Download
We download Uniprot search results for the parameters `(organism_id:9606) AND (reviewed:true)`, these can be accessed through a direct link: https://www.uniprot.org/uniprotkb?query=(organism_id:9606)%20AND%20(reviewed:true)

On the results page, click the *Download* button up top, then select *TSV* format from the drop-down menu. 
The default columns (`Entry`, `Reviewed`, `Entry Name`, `Protein names`, `Gene Names`, `Organism`) are fine, so all that's left is to download the file.

##### Data Preparation
The downloaded TSV file is ready to be used in the scripts, we saved it under the name `2022_03_uniprot_human_reviewed.tsv`, indicating the Uniprot release version `2022_03`, so that's the file name you will find in the scripts.

### Human Protein Atlas

##### Data Download
The *Human Protein Atlas* (HPA) project provides several datasets for direct download on their website ([Link](https://www.proteinatlas.org/about/download)).
We are interested in their TSV file with content corresponding to the data seen in HPA search results.
The most recent version can be downloaded directly from this url: https://www.proteinatlas.org/download/proteinatlas.tsv.zip

##### Data Preparation
After uncompressing the downloaded ZIP archive, the TSV file is ready to be used in the scripts. 
We have named the TSV file `v21.1-proteinatlas.tsv`, indicating the HPA release version `v21.1`.

### Pharos / TCRD
The bulk download option for the data accessible through the [Pharos](https://pharos.nih.gov/) platform is called the *Target Central Resource Database* (TCRD) and maintained by the University of New Mexico (UNM) School of Medicine [here](http://juniper.health.unm.edu/tcrd/).

##### Data Download
We can directly download a gz-compressed sql dump file from the UNM website.
In the first BALDR release, we used TCRD version 6.12.4 (released 2021-10-29), which can be downloaded via the url http://juniper.health.unm.edu/tcrd/download/old_versions/TCRDv6.12.4.sql.gz

While the file is only 7.5 GB in size, the download server is not the fastest, so downloading takes a bit of time. Once uncompressed, the sql dump file takes around 41 GB of storage space; loading this dump into an actually active database will take up another ca. 80 GB of storage.

##### Data Preparation
In order to get the TCRD data ready for R, we need to load the downloaded SQL dump into a local MySQL database which R can then connect to. I can recommend [this blog post](https://programminghistorian.org/en/lessons/getting-started-with-mysql-using-r) for familiarising yourself with the subject matter.

After installing and setting up a MySQL client, we can import the SQL dump file into a local database. First, we launch an interactive MySQL session in the terminal:
```
mysql --user=root --password=<your password> --host=localhost
```

In order to import the SQL dump file (a rather long task, although I managed to run it overnight on a Macbook Pro), run the following code inside the interactive MySQL session:
```
CREATE DATABASE tcrdimport;

USE tcrdimport;
SET autocommit=0 ; source <your/path/to/uncompressed_sql_file> ; COMMIT ;
```

The final preparation step that is required so that our R scripts can talk to the database is to set up a configuration file called `.my.cnf` in your home directory. 
Once you have set up this file, add an entry for the `tcrdimport` database that looks like this:
```
[tcrd]
database="tcrdimport"
user="root"
password="your password goes here"
```

### OpenTargets Overall JSON

##### Data Download
The *OpenTargets* platform provides various datasets for bulk download [here](https://platform.opentargets.org/downloads/data), along with a quick guide [here](https://platform-docs.opentargets.org/data-access/datasets). 

We use the UNIX utility `rsync` to download the *core annotation for targets* (version release 22.06, released 2022-06-24) data in JSON format:
```
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/22.06/output/etl/json/targets .
```

Note that this will download the set of target annotations, split into 200 smaller files, into a local directory (in this command the target directory is `.`, i.e. the current directory).

##### Data Preparation
In order to get the JSON files ready for our scripts, we concatenate the 200 downloaded files into one combined file called `combined_targets.json`, using the UNIX utility `cat`:
```
cat *.json > combined_targets.json
```

### OpenTargets Disease Associations
We will be downloading thirteen TSV files manually from the *OpenTargets* platform.

##### Data Download
We start off with the disease associations for the general *Diabetes mellitus* term which has the EFO ID `EFO_0000400`.  
We can find known Diabetes mellitus associations on this *OpenTargets* page: https://platform.opentargets.org/disease/EFO_0000400/associations

Here, we can download the displayed data by clicking on *Download table as TSV* in the top-right corner.
This starts the download of a file called `EFO_0000400-associated-diseases.tsv`.

We will repeat this procedure for the twelve selected diabetes subtypes in the following table:

|EFO ID          |Disease                                                       |
|:---------------|:-------------------------------------------------------------|
|Orphanet_183625 |Rare genetic diabetes mellitus                                |
|EFO_1001511     |monogenic diabetes                                            |
|EFO_1001121     |prediabetes syndrome                                          |
|EFO_0010164     |insulin-resistant diabetes mellitus                           |
|EFO_0007498     |Stiff-Person syndrome                                         |
|EFO_0004593     |gestational diabetes                                          |
|MONDO_0019193   |acquired generalized lipodystrophy                            |
|MONDO_0018575   |microcephalic primordial dwarfism-insulin resistance syndrome |
|MONDO_0014497   |polyendocrine-polyneuropathy syndrome                         |
|MONDO_0005148   |type 2 diabetes mellitus                                      |
|MONDO_0005147   |type 1 diabetes mellitus                                      |
|HP_0000857      |Neonatal insulin-dependent diabetes mellitus                  |

The URLs for disease association info of these twelve subtypes are as follows:

|URL                                                      |
|:--------------------------------------------------------|
|https://platform.opentargets.org/disease/Orphanet_183625/associations |
|https://platform.opentargets.org/disease/EFO_1001511/associations     |
|https://platform.opentargets.org/disease/EFO_1001121/associations     |
|https://platform.opentargets.org/disease/EFO_0010164/associations     |
|https://platform.opentargets.org/disease/EFO_0007498/associations     |
|https://platform.opentargets.org/disease/EFO_0004593/associations     |
|https://platform.opentargets.org/disease/MONDO_0019193/associations   |
|https://platform.opentargets.org/disease/MONDO_0018575/associations   |
|https://platform.opentargets.org/disease/MONDO_0014497/associations   |
|https://platform.opentargets.org/disease/MONDO_0005148/associations   |
|https://platform.opentargets.org/disease/MONDO_0005147/associations   |
|https://platform.opentargets.org/disease/HP_0000857/associations      |

##### Data Preparation
The thirteen TSV files are ready for use with our scripts without any changes.

### Drugbank

##### Data Download
DrugBank offers several non-commercial datasets for download on their website ([Link](https://go.drugbank.com/releases/latest)).
We are interested in their *Complete Database*, which we can access after a fairly quick registration process as academic user.

This complete dataset is available as a 145MB ZIP archive at the following URL: https://go.drugbank.com/releases/5-1-9/downloads/all-full-database

##### Data Preparation
After uncompressing the downloaded ZIP archive, we are presented with a 1.4GB XML file containing a lot of DrugBank's data. 

In order for this XML file to play nice with our scripts, there is one required preparation step: removing the namespace parameters in the beginning of the XML file.

The second line of the downloaded file contains a `<drugbank>` tag with a few extra parameters included. Overall, it looks like this:
```
<drugbank xmlns="http://www.drugbank.ca" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.drugbank.ca http://www.drugbank.ca/docs/drugbank.xsd" version="5.1" exported-on="2022-01-03">
```

We want to remove all the parameters from this tag, essentially changing the second line to only be:
```
<drugbank>
```

While this is a very straightforward change conceptually, of course it entails editing a fairly large file (and then saving the modified fairly large file), but with a modern computer and a bit of patience you should eventually be fine.

Once the `<drugbank>` tag is stripped of all parameters, the file is ready for use in our R scripts.
We have named it `2022-07-18-namespace_free_drugbank.xml`, indicating the download date (this corresponds to the version `5.1.9` released on 2022-01-04).

### Session Info
Our scripts load three libraries: `tidyverse`, `DBI`, and `xml2`.

```
> sessioninfo::session_info()
─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.0.5 (2021-03-31)
 os       macOS 12.6                  
 system   x86_64, darwin17.0          
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Copenhagen           
 date     2022-10-19                  

─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────
 package     * version date       lib source        
 assertthat    0.2.1   2019-03-21 [1] CRAN (R 4.0.0)
 backports     1.1.7   2020-05-13 [1] CRAN (R 4.0.0)
 bit           4.0.4   2020-08-04 [1] CRAN (R 4.0.2)
 bit64         4.0.5   2020-08-30 [1] CRAN (R 4.0.2)
 broom         0.7.6   2021-04-05 [1] CRAN (R 4.0.2)
 cellranger    1.1.0   2016-07-27 [1] CRAN (R 4.0.0)
 cli           2.5.0   2021-04-26 [1] CRAN (R 4.0.2)
 colorspace    2.0-1   2021-05-04 [1] CRAN (R 4.0.2)
 crayon        1.4.1   2021-02-08 [1] CRAN (R 4.0.2)
 DBI         * 1.1.0   2019-12-15 [1] CRAN (R 4.0.0)
 dbplyr        2.1.1   2021-04-06 [1] CRAN (R 4.0.2)
 digest        0.6.27  2020-10-24 [1] CRAN (R 4.0.2)
 dplyr       * 1.0.6   2021-05-05 [1] CRAN (R 4.0.2)
 ellipsis      0.3.2   2021-04-29 [1] CRAN (R 4.0.2)
 evaluate      0.16    2022-08-09 [1] CRAN (R 4.0.5)
 fansi         0.4.2   2021-01-15 [1] CRAN (R 4.0.2)
 fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.0.2)
 forcats     * 0.5.1   2021-01-27 [1] CRAN (R 4.0.2)
 fs            1.4.1   2020-04-04 [1] CRAN (R 4.0.0)
 generics      0.0.2   2018-11-29 [1] CRAN (R 4.0.0)
 ggplot2     * 3.3.5   2021-06-25 [1] CRAN (R 4.0.2)
 glue          1.4.2   2020-08-27 [1] CRAN (R 4.0.2)
 gtable        0.3.0   2019-03-25 [1] CRAN (R 4.0.0)
 haven         2.3.1   2020-06-01 [1] CRAN (R 4.0.0)
 hms           1.0.0   2021-01-13 [1] CRAN (R 4.0.2)
 htmltools     0.5.2   2021-08-25 [1] CRAN (R 4.0.2)
 httr          1.4.2   2020-07-20 [1] CRAN (R 4.0.2)
 jsonlite      1.7.2   2020-12-09 [1] CRAN (R 4.0.2)
 knitr         1.40    2022-08-24 [1] CRAN (R 4.0.5)
 lifecycle     1.0.0   2021-02-15 [1] CRAN (R 4.0.2)
 lubridate     1.7.10  2021-02-26 [1] CRAN (R 4.0.2)
 magrittr      2.0.1   2020-11-17 [1] CRAN (R 4.0.2)
 modelr        0.1.8   2020-05-19 [1] CRAN (R 4.0.0)
 munsell       0.5.0   2018-06-12 [1] CRAN (R 4.0.0)
 pillar        1.6.0   2021-04-13 [1] CRAN (R 4.0.0)
 pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.0.0)
 purrr       * 0.3.4   2020-04-17 [1] CRAN (R 4.0.0)
 R6            2.5.0   2020-10-28 [1] CRAN (R 4.0.2)
 Rcpp          1.0.7   2021-07-07 [1] CRAN (R 4.0.2)
 readr       * 1.4.0   2020-10-05 [1] CRAN (R 4.0.2)
 readxl        1.3.1   2019-03-13 [1] CRAN (R 4.0.0)
 reprex        2.0.0   2021-04-02 [1] CRAN (R 4.0.2)
 rlang         0.4.11  2021-04-30 [1] CRAN (R 4.0.2)
 RMariaDB      1.2.1   2021-12-20 [1] CRAN (R 4.0.2)
 rmarkdown     2.16    2022-08-24 [1] CRAN (R 4.0.5)
 rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.0.2)
 rvest         1.0.0   2021-03-09 [1] CRAN (R 4.0.2)
 scales        1.1.1   2020-05-11 [1] CRAN (R 4.0.0)
 sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 4.0.0)
 stringi       1.4.6   2020-02-17 [1] CRAN (R 4.0.0)
 stringr     * 1.4.0   2019-02-10 [1] CRAN (R 4.0.0)
 tibble      * 3.1.1   2021-04-18 [1] CRAN (R 4.0.2)
 tidyr       * 1.1.3   2021-03-03 [1] CRAN (R 4.0.2)
 tidyselect    1.1.0   2020-05-11 [1] CRAN (R 4.0.0)
 tidyverse   * 1.3.1   2021-04-15 [1] CRAN (R 4.0.2)
 utf8          1.2.1   2021-03-12 [1] CRAN (R 4.0.2)
 vctrs         0.3.8   2021-04-29 [1] CRAN (R 4.0.2)
 withr         2.4.2   2021-04-18 [1] CRAN (R 4.0.2)
 xfun          0.30    2022-03-02 [1] CRAN (R 4.0.5)
 xml2        * 1.3.3   2021-11-30 [1] CRAN (R 4.0.2)
```
