---
author : "L. Cauquil"
title: "Modèle linéaire mixte avec les packages purrr et broom.mixed"
date: "26-01-2022"
output:
  html_document: 
    code_folding: show
    toc: yes
    toc_float: yes
    keep_md: TRUE
    css: style.css
  pdf_document:
    toc: yes
  word_document:
    toc: yes
editor_options: 
  chunk_output_type: inline
---



## **Packages et importation des données**

### **Packages**


```r
## Manipulations des données
library(tibble)
library(dplyr)
library(tidyr)

## Mise en forme des résultats stat
library(purrr)
library(broom)
library(broom.mixed)

## Objet phyloseq
library(phyloseq)

## Fonctions stat
library(car)
library(multcomp)
library(lmerTest)
library(emmeans)

## Présentation des tables
library(DT)
library(flextable)
# library(xlsx)
```

### **Importation des données**

On utilise l'objet phyloseq `tab_phylum`


```r
load("data/data_phylum.RData")
tab_Phylum
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 9 taxa and 40 samples ]
## sample_data() Sample Data:       [ 40 samples by 9 sample variables ]
## tax_table()   Taxonomy Table:    [ 9 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 9 tips and 8 internal nodes ]
```

Il y a 9 phyla présents dans l'écosystème.


## **Construction des objets**

Différents objets créés:

 - `RecapPhylum`: data.frame avec les moyennes, sd et SEM des abondances relatives de chaque phylum sans tenir compte des groupes
 - `top_Phylum`: vecteur de Phyla qui ont au moins une moyenne des abondances relatives au sein d'un groupe Segment_Age > à 0,005 
 - Rajout de la variable `Firm_bact_ratio` = ration Firmicutes/Bacteroidota
 - `df_stat`: data.frame utilisé pour appliqué un modèle linéaire mixte sur les abondances relatives des phyla sélectionnés + `Firm_bact_ratio`

Les analyses sont faites sur les abondances relatives


```r
## Transforme en abondance relative
relatabunphy <- transform_sample_counts(tab_Phylum, function(OTU) OTU/sum(OTU))
```

La fonction `psmelt()` fusionne les tables `otu_table`, `sample_otu` et `tax_table` d'un objet phyloseq et crée un seul data.frame au format long.  
Les phyla sont regroupés dans une seule colonne Phylum !


```r
dfrap <- psmelt(relatabunphy) 
head(dfrap)
```

```
##          OTU                      Sample Abundance       GPS Piglet Age Sexe
## 39 Cluster_1 GPS160115_ATGAAC-JLGMN_L001 0.9995071 GPS160115     14 D35    F
## 38 Cluster_1 GPS161671_CTCTAC-JLGMN_L001 0.9993956 GPS161671     23 D35    M
## 19 Cluster_1 GPS158690_CCTTGA-JLGMN_L001 0.9984081 GPS158690      5 D35    M
## 22 Cluster_1 GPS158965_CACCCA-JLGMN_L001 0.9983199 GPS158965      9 D21    F
## 18 Cluster_1 GPS161675_ACCGTG-JLGMN_L001 0.9981725 GPS161675     24 D35    M
## 26 Cluster_1 GPS160464_CCCAAA-JLGMN_L001 0.9978882 GPS160464     20 D21    F
##    Segment Segment_Age Concentration.ADN                    SampleID
## 39 Jejunum Jejunum_D35             7,736 GPS160115_ATGAAC-JLGMN_L001
## 38 Jejunum Jejunum_D35              9,75 GPS161671_CTCTAC-JLGMN_L001
## 19 Jejunum Jejunum_D35             75,69 GPS158690_CCTTGA-JLGMN_L001
## 22 Jejunum Jejunum_D21             25,34 GPS158965_CACCCA-JLGMN_L001
## 18 Jejunum Jejunum_D35             5,389 GPS161675_ACCGTG-JLGMN_L001
## 26 Jejunum Jejunum_D21             18,19 GPS160464_CCCAAA-JLGMN_L001
##    Slaughter_date  Kingdom     Phylum
## 39             D3 Bacteria Firmicutes
## 38             D5 Bacteria Firmicutes
## 19             D1 Bacteria Firmicutes
## 22             D2 Bacteria Firmicutes
## 18             D5 Bacteria Firmicutes
## 26             D4 Bacteria Firmicutes
```

### **Tableau général avec mean, sd, SEM des abondances relatives par phylum**

Construction de la table: tableau général sans tenir des groupes  


```r
RecapPhylum <- dfrap |>  
  dplyr::select(Phylum,Abundance) |> 
  group_by(Phylum) |> 
  summarise(data.frame(mean = mean(Abundance),
                       sd = sd(Abundance), 
                       sem = sd(Abundance) / sqrt(length(Abundance))))
RecapPhylum[,2:4] <- round(RecapPhylum[,2:4],4)*100

datatable(RecapPhylum)
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-b4d8925adc6cf284c63e" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-b4d8925adc6cf284c63e">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9"],["Actinobacteriota","Bacteroidota","Campylobacterota","Desulfobacterota","Firmicutes","Fusobacteriota","Patescibacteria","Proteobacteria","Spirochaetota"],[0.43,11.58,0.11,0.55,82.87,0.81,0.06,3.56,0.03],[0.79,14.56,0.27,0.88,18.7,2.55,0.3,8.27,0.09],[0.12,2.3,0.04,0.14,2.96,0.4,0.05,1.31,0.01]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Phylum<\/th>\n      <th>mean<\/th>\n      <th>sd<\/th>\n      <th>sem<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


### **Sélection des phyla avec au moins une moyenne d'abondance relative > à 0,005 au sein d'un groupe Segment_Age**
  
Sélection des phyla  


```r
top_Phylum <- dfrap |>  
  dplyr::select(Phylum, Abundance, Segment_Age) |> 
  group_by(Phylum,Segment_Age) |> 
  summarise_all(list(mean = mean)) |> 
  dplyr::arrange(desc(mean)) |> 
  dplyr::filter(mean >= 0.005) |> 
  dplyr::arrange(Phylum, Segment_Age) |> 
  distinct(Phylum) |> 
  as_vector()
top_Phylum
```

```
##            Phylum1            Phylum2            Phylum3            Phylum4 
## "Actinobacteriota"     "Bacteroidota" "Desulfobacterota"       "Firmicutes" 
##            Phylum5            Phylum6 
##   "Fusobacteriota"   "Proteobacteria"
```

Selection des phyla dans la table générale


```r
dfrap <- dfrap |> 
  dplyr::filter(Phylum %in% top_Phylum)
```

### **Rajout de la variable `Firm_bact_ratio` aux phyla**


```r
dfrap |> 
  dplyr::select(Phylum, Abundance, Segment_Age, Piglet, Slaughter_date, Segment, Age) |> 
#  arrange(Phylum) |> 
  pivot_wider(names_from = Phylum,
              values_from = Abundance) |> 
  mutate(Firm_bact_ratio = Firmicutes/Bacteroidota) |> 
  pivot_longer(cols = where(is.numeric),
               names_to = "phylum", 
               values_to = "abundance") -> df_stat
```

Transformation des abundances en racine de racine


```r
df_stat$abundance <- df_stat$abundance^0.25
```

## **Statistiques**  

**Paramètres des contrastes**


```r
options(contrasts = c("contr.sum", "contr.poly"))
```

**Modèle linéaire mixte**

Le modèle utilise 

- la fonction `lmer` du package `lmerTest` 
- l'effet piglet en variable aléatoire (1|Piglet) 

Mais Piglet est niché dans slaughter_date

- l'addition de Slaughter_date sur (1|Slaughter_date/Piglet)

NB l'écriture (1|Slaughter_date/Piglet) est équivalente à (1|Piglet) +(1|Slaughter_date) car les porcs ont des identifiants différents entre les date d'abattage
plus d'explication ici https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified  

https://github.com/lcauquil/tuto_mix_model/blob/master/Modele_mixte_guidelines_2019.pdf pour les explications REML (REML pour restricted maximum likelihood) et ML (maximum likelihood)

### **Modèle utilisé**

**lmer(abundance ~ Segment * Age + (1|Piglet)+ (1|Slaughter_date), data = as.data.frame(data))**

Application du modèle sur le phylum "Bacteroidota"


```r
df_stat |> 
  filter(phylum == "Bacteroidota") -> lm_bacte

lmer(abundance ~ Segment * Age + (1|Piglet)+ (1|Slaughter_date), data = lm_bacte)
```

```
## Linear mixed model fit by REML ['lmerModLmerTest']
## Formula: abundance ~ Segment * Age + (1 | Piglet) + (1 | Slaughter_date)
##    Data: lm_bacte
## REML criterion at convergence: -31.5725
## Random effects:
##  Groups         Name        Std.Dev.
##  Piglet         (Intercept) 0.14663 
##  Slaughter_date (Intercept) 0.00000 
##  Residual                   0.06558 
## Number of obs: 40, groups:  Piglet, 24; Slaughter_date, 6
## Fixed Effects:
##   (Intercept)       Segment1           Age1  Segment1:Age1  
##      0.463270      -0.266397       0.026846      -0.004721  
## optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings
```

### **Fonctions**

On crées plusieurs fonctions à appliquer pour chacun des phylum.  
Fonctions pour :

 - appliquer le modèle
 - effectuer le test shapiro
 - calcluer les p-values
 - récupérer les lettres

Pour appliquer chaque fonction on utilise la fonction `map()` du package purrr. Les fonctions `map()` et dérivées visent à remplacer les boucles et les fonctions de la famille des apply.  

Seule contrainte des fonctions `map()`, elles s'appliquent uniquement à des listes. Il faut donc convertir les data.frames en liste, mais c'est facile avec la fonction `split()`.

**Exemple d'utilisation de la fonction `map()`**

La fonction `map()` à généralement 2 arguments: une liste et une fonction


```r
x <- list(a = 1:10, b = 5:110)
x
```

```
## $a
##  [1]  1  2  3  4  5  6  7  8  9 10
## 
## $b
##   [1]   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22
##  [19]  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40
##  [37]  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58
##  [55]  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76
##  [73]  77  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94
##  [91]  95  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110
```

```r
map(x, mean)
```

```
## $a
## [1] 5.5
## 
## $b
## [1] 57.5
```

```r
## sortie sous forme de data.frame
map_df(x, mean)
```

```
## # A tibble: 1 × 2
##       a     b
##   <dbl> <dbl>
## 1   5.5  57.5
```

```r
rm(x)
```

**Fonction split()**

La fonction `split()` sépare un data.frame en liste à partir des différents niveaux d'une variable
 

```r
split(df_stat, df_stat$phylum)
```

```{.scroll-300}
## $Actinobacteriota
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum           abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>                <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Actinobacteriota     0.115
##  2 Jejunum_D35 23     D5             Jejunum D35   Actinobacteriota     0.109
##  3 Jejunum_D35 5      D1             Jejunum D35   Actinobacteriota     0.142
##  4 Jejunum_D21 9      D2             Jejunum D21   Actinobacteriota     0.136
##  5 Jejunum_D35 24     D5             Jejunum D35   Actinobacteriota     0.145
##  6 Jejunum_D21 20     D4             Jejunum D21   Actinobacteriota     0.120
##  7 Jejunum_D21 27     D6             Jejunum D21   Actinobacteriota     0.203
##  8 Jejunum_D21 26     D6             Jejunum D21   Actinobacteriota     0.178
##  9 Jejunum_D21 18     D4             Jejunum D21   Actinobacteriota     0.149
## 10 Jejunum_D35 7      D1             Jejunum D35   Actinobacteriota     0.125
## # … with 30 more rows
## 
## $Bacteroidota
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum       abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>            <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Bacteroidota    0.0966
##  2 Jejunum_D35 23     D5             Jejunum D35   Bacteroidota    0.120 
##  3 Jejunum_D35 5      D1             Jejunum D35   Bacteroidota    0.129 
##  4 Jejunum_D21 9      D2             Jejunum D21   Bacteroidota    0.139 
##  5 Jejunum_D35 24     D5             Jejunum D35   Bacteroidota    0.150 
##  6 Jejunum_D21 20     D4             Jejunum D21   Bacteroidota    0.131 
##  7 Jejunum_D21 27     D6             Jejunum D21   Bacteroidota    0.0855
##  8 Jejunum_D21 26     D6             Jejunum D21   Bacteroidota    0.141 
##  9 Jejunum_D21 18     D4             Jejunum D21   Bacteroidota    0.120 
## 10 Jejunum_D35 7      D1             Jejunum D35   Bacteroidota    0.121 
## # … with 30 more rows
## 
## $Desulfobacterota
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum           abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>                <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Desulfobacterota         0
##  2 Jejunum_D35 23     D5             Jejunum D35   Desulfobacterota         0
##  3 Jejunum_D35 5      D1             Jejunum D35   Desulfobacterota         0
##  4 Jejunum_D21 9      D2             Jejunum D21   Desulfobacterota         0
##  5 Jejunum_D35 24     D5             Jejunum D35   Desulfobacterota         0
##  6 Jejunum_D21 20     D4             Jejunum D21   Desulfobacterota         0
##  7 Jejunum_D21 27     D6             Jejunum D21   Desulfobacterota         0
##  8 Jejunum_D21 26     D6             Jejunum D21   Desulfobacterota         0
##  9 Jejunum_D21 18     D4             Jejunum D21   Desulfobacterota         0
## 10 Jejunum_D35 7      D1             Jejunum D35   Desulfobacterota         0
## # … with 30 more rows
## 
## $Firm_bact_ratio
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum          abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>               <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Firm_bact_ratio     10.4 
##  2 Jejunum_D35 23     D5             Jejunum D35   Firm_bact_ratio      8.31
##  3 Jejunum_D35 5      D1             Jejunum D35   Firm_bact_ratio      7.73
##  4 Jejunum_D21 9      D2             Jejunum D21   Firm_bact_ratio      7.20
##  5 Jejunum_D35 24     D5             Jejunum D35   Firm_bact_ratio      6.65
##  6 Jejunum_D21 20     D4             Jejunum D21   Firm_bact_ratio      7.63
##  7 Jejunum_D21 27     D6             Jejunum D21   Firm_bact_ratio     11.7 
##  8 Jejunum_D21 26     D6             Jejunum D21   Firm_bact_ratio      7.10
##  9 Jejunum_D21 18     D4             Jejunum D21   Firm_bact_ratio      8.35
## 10 Jejunum_D35 7      D1             Jejunum D35   Firm_bact_ratio      8.29
## # … with 30 more rows
## 
## $Firmicutes
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum     abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>          <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Firmicutes     1.00 
##  2 Jejunum_D35 23     D5             Jejunum D35   Firmicutes     1.00 
##  3 Jejunum_D35 5      D1             Jejunum D35   Firmicutes     1.00 
##  4 Jejunum_D21 9      D2             Jejunum D21   Firmicutes     1.00 
##  5 Jejunum_D35 24     D5             Jejunum D35   Firmicutes     1.00 
##  6 Jejunum_D21 20     D4             Jejunum D21   Firmicutes     0.999
##  7 Jejunum_D21 27     D6             Jejunum D21   Firmicutes     0.999
##  8 Jejunum_D21 26     D6             Jejunum D21   Firmicutes     0.999
##  9 Jejunum_D21 18     D4             Jejunum D21   Firmicutes     0.999
## 10 Jejunum_D35 7      D1             Jejunum D35   Firmicutes     0.999
## # … with 30 more rows
## 
## $Fusobacteriota
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum         abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>              <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Fusobacteriota    0     
##  2 Jejunum_D35 23     D5             Jejunum D35   Fusobacteriota    0     
##  3 Jejunum_D35 5      D1             Jejunum D35   Fusobacteriota    0     
##  4 Jejunum_D21 9      D2             Jejunum D21   Fusobacteriota    0     
##  5 Jejunum_D35 24     D5             Jejunum D35   Fusobacteriota    0     
##  6 Jejunum_D21 20     D4             Jejunum D21   Fusobacteriota    0.120 
##  7 Jejunum_D21 27     D6             Jejunum D21   Fusobacteriota    0.0855
##  8 Jejunum_D21 26     D6             Jejunum D21   Fusobacteriota    0     
##  9 Jejunum_D21 18     D4             Jejunum D21   Fusobacteriota    0.0846
## 10 Jejunum_D35 7      D1             Jejunum D35   Fusobacteriota    0.0770
## # … with 30 more rows
## 
## $Proteobacteria
## # A tibble: 40 × 7
##    Segment_Age Piglet Slaughter_date Segment Age   phylum         abundance
##    <fct>       <fct>  <chr>          <fct>   <fct> <chr>              <dbl>
##  1 Jejunum_D35 14     D3             Jejunum D35   Proteobacteria     0.123
##  2 Jejunum_D35 23     D5             Jejunum D35   Proteobacteria     0.126
##  3 Jejunum_D35 5      D1             Jejunum D35   Proteobacteria     0.172
##  4 Jejunum_D21 9      D2             Jejunum D21   Proteobacteria     0.176
##  5 Jejunum_D35 24     D5             Jejunum D35   Proteobacteria     0.172
##  6 Jejunum_D21 20     D4             Jejunum D21   Proteobacteria     0.194
##  7 Jejunum_D21 27     D6             Jejunum D21   Proteobacteria     0.141
##  8 Jejunum_D21 26     D6             Jejunum D21   Proteobacteria     0.167
##  9 Jejunum_D21 18     D4             Jejunum D21   Proteobacteria     0.199
## 10 Jejunum_D35 7      D1             Jejunum D35   Proteobacteria     0.218
## # … with 30 more rows
```

**Fonction du modèle**


```r
fit_model <- function(data) {
  mod_REML <- lmer(abundance ~ Segment * Age + (1|Piglet)+ (1|Slaughter_date), data = as.data.frame(data))
	mod_REML_update <- update(mod_REML, REML = F)
}
```

**Application de la fonction `fit_model` sur les Firmicutes**  


```r
df_stat |> 
  dplyr::filter(phylum == "Firmicutes") |> 
  fit_model() |> 
  summary()
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's
##   method [lmerModLmerTest]
## Formula: abundance ~ Segment * Age + (1 | Piglet) + (1 | Slaughter_date)
##    Data: as.data.frame(data)
## 
##      AIC      BIC   logLik deviance df.resid 
##   -108.8    -97.0     61.4   -122.8       33 
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -1.79409 -0.24396  0.06952  0.20530  1.87531 
## 
## Random effects:
##  Groups         Name        Variance  Std.Dev.
##  Piglet         (Intercept) 0.0045381 0.06737 
##  Slaughter_date (Intercept) 0.0000000 0.00000 
##  Residual                   0.0005582 0.02363 
## Number of obs: 40, groups:  Piglet, 24; Slaughter_date, 6
## 
## Fixed effects:
##                Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)    0.935449   0.014454 23.079580  64.720  < 2e-16 ***
## Segment1       0.034004   0.004453 14.911031   7.637 1.58e-06 ***
## Age1          -0.016489   0.014454 23.079580  -1.141    0.266    
## Segment1:Age1  0.004979   0.004453 14.911031   1.118    0.281    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) Sgmnt1 Age1  
## Segment1    -0.127              
## Age1        -0.034  0.112       
## Segmnt1:Ag1  0.112 -0.363 -0.127
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

**Application de la fonction `fit_model` sur tous les phyla**  


```r
df_stat |> 
  split(df_stat$phylum) |> 
  map(fit_model) -> result
```

**Resultat pour les Firmicutes**


```r
result$Firmicutes |> summary()
```

```
## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's
##   method [lmerModLmerTest]
## Formula: abundance ~ Segment * Age + (1 | Piglet) + (1 | Slaughter_date)
##    Data: as.data.frame(data)
## 
##      AIC      BIC   logLik deviance df.resid 
##   -108.8    -97.0     61.4   -122.8       33 
## 
## Scaled residuals: 
##      Min       1Q   Median       3Q      Max 
## -1.79409 -0.24396  0.06952  0.20530  1.87531 
## 
## Random effects:
##  Groups         Name        Variance  Std.Dev.
##  Piglet         (Intercept) 0.0045381 0.06737 
##  Slaughter_date (Intercept) 0.0000000 0.00000 
##  Residual                   0.0005582 0.02363 
## Number of obs: 40, groups:  Piglet, 24; Slaughter_date, 6
## 
## Fixed effects:
##                Estimate Std. Error        df t value Pr(>|t|)    
## (Intercept)    0.935449   0.014454 23.079580  64.720  < 2e-16 ***
## Segment1       0.034004   0.004453 14.911031   7.637 1.58e-06 ***
## Age1          -0.016489   0.014454 23.079580  -1.141    0.266    
## Segment1:Age1  0.004979   0.004453 14.911031   1.118    0.281    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Correlation of Fixed Effects:
##             (Intr) Sgmnt1 Age1  
## Segment1    -0.127              
## Age1        -0.034  0.112       
## Segmnt1:Ag1  0.112 -0.363 -0.127
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

### **Mise en forme des résultats**

Pour mettre en forme les résultats on peut utiliser la fonction `tidy()` des packages broom et broom.mixed qui récupèrent les résultats des tests et les formatent en data.frame, toujours en utilisant la fonction `map()`

La fonction `tidy()` s'applique à de nombreux tests:

[Listes des méthodes prises en charge par broom](https://broom.tidymodels.org/articles/available-methods.html){target="blank"}

[Listes des méthodes prises en charge par broom.mixed](https://cran.r-project.org/web/packages/broom.mixed/vignettes/broom_mixed_intro.html){target="blank"}

Mise en forme des résultats pour les Firmicutes


```r
tidy(result$Firmicutes)
```

```
## # A tibble: 7 × 8
##   effect   group          term         estimate std.er…¹ stati…²    df   p.value
##   <chr>    <chr>          <chr>           <dbl>    <dbl>   <dbl> <dbl>     <dbl>
## 1 fixed    <NA>           (Intercept)   0.935    0.0145    64.7   23.1  1.33e-27
## 2 fixed    <NA>           Segment1      0.0340   0.00445    7.64  14.9  1.58e- 6
## 3 fixed    <NA>           Age1         -0.0165   0.0145    -1.14  23.1  2.66e- 1
## 4 fixed    <NA>           Segment1:Ag…  0.00498  0.00445    1.12  14.9  2.81e- 1
## 5 ran_pars Piglet         sd__(Interc…  0.0674  NA         NA     NA   NA       
## 6 ran_pars Slaughter_date sd__(Interc…  0       NA         NA     NA   NA       
## 7 ran_pars Residual       sd__Observa…  0.0236  NA         NA     NA   NA       
## # … with abbreviated variable names ¹​std.error, ²​statistic
```

Appliqué à tous les phyla


```r
result_model <- df_stat |>
  split(df_stat$phylum) |> 
  map(fit_model) |> 
  map(broom.mixed::tidy, .id = "phylum")
result_model
```

```{.scroll-300}
## $Actinobacteriota
## # A tibble: 7 × 8
##   effect   group          term         estimate std.er…¹ stati…²    df   p.value
##   <chr>    <chr>          <chr>           <dbl>    <dbl>   <dbl> <dbl>     <dbl>
## 1 fixed    <NA>           (Intercept)   0.220    0.0145   15.1    21.7  5.34e-13
## 2 fixed    <NA>           Segment1     -0.0124   0.00974  -1.27   15.3  2.23e- 1
## 3 fixed    <NA>           Age1          0.0130   0.0145    0.896  21.7  3.80e- 1
## 4 fixed    <NA>           Segment1:Ag…  0.00247  0.00974   0.254  15.3  8.03e- 1
## 5 ran_pars Piglet         sd__(Interc…  0.0527  NA        NA      NA   NA       
## 6 ran_pars Slaughter_date sd__(Interc…  0       NA        NA      NA   NA       
## 7 ran_pars Residual       sd__Observa…  0.0541  NA        NA      NA   NA       
## # … with abbreviated variable names ¹​std.error, ²​statistic
## 
## $Bacteroidota
## # A tibble: 7 × 8
##   effect   group          term          estimate std.e…¹ stati…²    df   p.value
##   <chr>    <chr>          <chr>            <dbl>   <dbl>   <dbl> <dbl>     <dbl>
## 1 fixed    <NA>           (Intercept)    0.463    0.0310  15.0    22.6  3.37e-13
## 2 fixed    <NA>           Segment1      -0.267    0.0115 -23.2    14.4  7.60e-13
## 3 fixed    <NA>           Age1           0.0269   0.0310   0.867  22.6  3.95e- 1
## 4 fixed    <NA>           Segment1:Age1 -0.00476  0.0115  -0.414  14.4  6.85e- 1
## 5 ran_pars Piglet         sd__(Interce…  0.141   NA       NA      NA   NA       
## 6 ran_pars Slaughter_date sd__(Interce…  0       NA       NA      NA   NA       
## 7 ran_pars Residual       sd__Observat…  0.0612  NA       NA      NA   NA       
## # … with abbreviated variable names ¹​std.error, ²​statistic
## 
## $Desulfobacterota
## # A tibble: 7 × 8
##   effect   group          term           estimate std.e…¹ stati…²    df  p.value
##   <chr>    <chr>          <chr>             <dbl>   <dbl>   <dbl> <dbl>    <dbl>
## 1 fixed    <NA>           (Intercept)     0.180    0.0114  15.8    11.8  2.55e-9
## 2 fixed    <NA>           Segment1       -0.141    0.0109 -12.9    11.1  4.64e-8
## 3 fixed    <NA>           Age1            0.00305  0.0114   0.269  11.8  7.93e-1
## 4 fixed    <NA>           Segment1:Age1  -0.00312  0.0109  -0.286  11.1  7.80e-1
## 5 ran_pars Piglet         sd__(Intercep…  0.0157  NA       NA      NA   NA      
## 6 ran_pars Slaughter_date sd__(Intercep…  0       NA       NA      NA   NA      
## 7 ran_pars Residual       sd__Observati…  0.0640  NA       NA      NA   NA      
## # … with abbreviated variable names ¹​std.error, ²​statistic
## 
## $Firm_bact_ratio
## # A tibble: 7 × 8
##   effect   group          term            estim…¹ std.e…² stati…³    df  p.value
##   <chr>    <chr>          <chr>             <dbl>   <dbl>   <dbl> <dbl>    <dbl>
## 1 fixed    <NA>           (Intercept)       4.00    0.390  10.2   11.6   3.71e-7
## 2 fixed    <NA>           Segment1          2.79    0.289   9.64   7.75  1.39e-5
## 3 fixed    <NA>           Age1             -0.226   0.390  -0.579 11.6   5.74e-1
## 4 fixed    <NA>           Segment1:Age1    -0.106   0.289  -0.365  7.75  7.25e-1
## 5 ran_pars Piglet         sd__(Intercept)   1.28   NA      NA     NA    NA      
## 6 ran_pars Slaughter_date sd__(Intercept)   0      NA      NA     NA    NA      
## 7 ran_pars Residual       sd__Observation   1.63   NA      NA     NA    NA      
## # … with abbreviated variable names ¹​estimate, ²​std.error, ³​statistic
## 
## $Firmicutes
## # A tibble: 7 × 8
##   effect   group          term         estimate std.er…¹ stati…²    df   p.value
##   <chr>    <chr>          <chr>           <dbl>    <dbl>   <dbl> <dbl>     <dbl>
## 1 fixed    <NA>           (Intercept)   0.935    0.0145    64.7   23.1  1.33e-27
## 2 fixed    <NA>           Segment1      0.0340   0.00445    7.64  14.9  1.58e- 6
## 3 fixed    <NA>           Age1         -0.0165   0.0145    -1.14  23.1  2.66e- 1
## 4 fixed    <NA>           Segment1:Ag…  0.00498  0.00445    1.12  14.9  2.81e- 1
## 5 ran_pars Piglet         sd__(Interc…  0.0674  NA         NA     NA   NA       
## 6 ran_pars Slaughter_date sd__(Interc…  0       NA         NA     NA   NA       
## 7 ran_pars Residual       sd__Observa…  0.0236  NA         NA     NA   NA       
## # … with abbreviated variable names ¹​std.error, ²​statistic
## 
## $Fusobacteriota
## # A tibble: 7 × 8
##   effect   group          term            estim…¹ std.e…² stati…³    df  p.value
##   <chr>    <chr>          <chr>             <dbl>   <dbl>   <dbl> <dbl>    <dbl>
## 1 fixed    <NA>           (Intercept)      0.132   0.0272    4.85  22.5  7.18e-5
## 2 fixed    <NA>           Segment1        -0.0262  0.0157   -1.67  15.1  1.15e-1
## 3 fixed    <NA>           Age1             0.0854  0.0272    3.14  22.5  4.68e-3
## 4 fixed    <NA>           Segment1:Age1   -0.0339  0.0157   -2.17  15.1  4.65e-2
## 5 ran_pars Piglet         sd__(Intercept)  0.109  NA        NA     NA   NA      
## 6 ran_pars Slaughter_date sd__(Intercept)  0      NA        NA     NA   NA      
## 7 ran_pars Residual       sd__Observation  0.0857 NA        NA     NA   NA      
## # … with abbreviated variable names ¹​estimate, ²​std.error, ³​statistic
## 
## $Proteobacteria
## # A tibble: 7 × 8
##   effect   group          term            estim…¹ std.e…² stati…³    df  p.value
##   <chr>    <chr>          <chr>             <dbl>   <dbl>   <dbl> <dbl>    <dbl>
## 1 fixed    <NA>           (Intercept)     3.15e-1  0.0366  8.61    5.37  2.41e-4
## 2 fixed    <NA>           Segment1        4.97e-4  0.0211  0.0236 12.2   9.82e-1
## 3 fixed    <NA>           Age1            1.65e-2  0.0366  0.451   5.37  6.69e-1
## 4 fixed    <NA>           Segment1:Age1   5.62e-4  0.0211  0.0266 12.2   9.79e-1
## 5 ran_pars Piglet         sd__(Intercept) 1.23e-1 NA      NA      NA    NA      
## 6 ran_pars Slaughter_date sd__(Intercept) 3.93e-2 NA      NA      NA    NA      
## 7 ran_pars Residual       sd__Observation 1.16e-1 NA      NA      NA    NA      
## # … with abbreviated variable names ¹​estimate, ²​std.error, ³​statistic
```

Pour regrouper les résultats de chaque phylum dans une même table on utilise la function `map_df()` qui crée un directement un data.frame.


```r
result_model <- df_stat |>
  split(df_stat$phylum) |> 
  map(fit_model) |> 
  map_df(broom.mixed::tidy, .id = "phylum")
result_model
```

```
## # A tibble: 49 × 9
##    phylum           effect group term  estimate std.er…¹ stati…²    df   p.value
##    <chr>            <chr>  <chr> <chr>    <dbl>    <dbl>   <dbl> <dbl>     <dbl>
##  1 Actinobacteriota fixed  <NA>  (Int…  0.220    0.0145   15.1    21.7  5.34e-13
##  2 Actinobacteriota fixed  <NA>  Segm… -0.0124   0.00974  -1.27   15.3  2.23e- 1
##  3 Actinobacteriota fixed  <NA>  Age1   0.0130   0.0145    0.896  21.7  3.80e- 1
##  4 Actinobacteriota fixed  <NA>  Segm…  0.00247  0.00974   0.254  15.3  8.03e- 1
##  5 Actinobacteriota ran_p… Pigl… sd__…  0.0527  NA        NA      NA   NA       
##  6 Actinobacteriota ran_p… Slau… sd__…  0       NA        NA      NA   NA       
##  7 Actinobacteriota ran_p… Resi… sd__…  0.0541  NA        NA      NA   NA       
##  8 Bacteroidota     fixed  <NA>  (Int…  0.463    0.0310   15.0    22.6  3.37e-13
##  9 Bacteroidota     fixed  <NA>  Segm… -0.267    0.0115  -23.2    14.4  7.60e-13
## 10 Bacteroidota     fixed  <NA>  Age1   0.0269   0.0310    0.867  22.6  3.95e- 1
## # … with 39 more rows, and abbreviated variable names ¹​std.error, ²​statistic
```

### **Même méthode pour les p.value, shapiro et lettres**

**p.value**


```r
## fonction
p_val <- function(data)
{
  Anova(data, type = "III")
}

## résultats
result_pval <- df_stat |>
  split(df_stat$phylum) |> 
  map(fit_model) |>
  map(p_val) |> 
  map_dfr(tidy, .id = "phylum")
result_pval
```

```
## # A tibble: 28 × 5
##    phylum           term        statistic    df   p.value
##    <chr>            <chr>           <dbl> <dbl>     <dbl>
##  1 Actinobacteriota (Intercept)  229.         1 8.68e- 52
##  2 Actinobacteriota Segment        1.62       1 2.04e-  1
##  3 Actinobacteriota Age            0.803      1 3.70e-  1
##  4 Actinobacteriota Segment:Age    0.0645     1 7.99e-  1
##  5 Bacteroidota     (Intercept)  224.         1 1.52e- 50
##  6 Bacteroidota     Segment      540.         1 2.23e-119
##  7 Bacteroidota     Age            0.752      1 3.86e-  1
##  8 Bacteroidota     Segment:Age    0.172      1 6.79e-  1
##  9 Desulfobacterota (Intercept)  250.         1 2.48e- 56
## 10 Desulfobacterota Segment      168.         1 2.38e- 38
## # … with 18 more rows
```

**Ajout des p-adjusted**


```r
result_pval |> 
  #dplyr::filter(effect == "fixed") |> 
  dplyr::select(phylum, term, p.value) |> 
  pivot_wider(names_from = term,
              values_from = p.value) |> 
  mutate(p_adj_Segment = p.adjust(Segment, method = "BH"),
         p_adj_Age = p.adjust(Age, method = "BH"),
         p_adj_Age_Segment = p.adjust(`Segment:Age`, method = "BH")) -> df_pval
df_pval
```

```
## # A tibble: 7 × 8
##   phylum           (Interc…¹   Segment     Age Segme…² p_adj_S…³ p_adj…⁴ p_adj…⁵
##   <chr>                <dbl>     <dbl>   <dbl>   <dbl>     <dbl>   <dbl>   <dbl>
## 1 Actinobacteriota  8.68e-52 2.04e-  1 0.370    0.799  2.38e-  1  0.675    0.933
## 2 Bacteroidota      1.52e-50 2.23e-119 0.386    0.679  1.56e-118  0.675    0.933
## 3 Desulfobacterota  2.48e-56 2.38e- 38 0.788    0.775  8.33e- 38  0.788    0.933
## 4 Firm_bact_ratio   1.31e-24 5.42e- 22 0.563    0.715  1.27e- 21  0.760    0.933
## 5 Firmicutes        0        2.22e- 14 0.254    0.264  3.89e- 14  0.675    0.922
## 6 Fusobacteriota    1.24e- 6 9.45e-  2 0.00169  0.0302 1.32e-  1  0.0118   0.211
## 7 Proteobacteria    7.01e-18 9.81e-  1 0.652    0.979  9.81e-  1  0.760    0.979
## # … with abbreviated variable names ¹​`(Intercept)`, ²​`Segment:Age`,
## #   ³​p_adj_Segment, ⁴​p_adj_Age, ⁵​p_adj_Age_Segment
```

**Shapiro**


```r
## fonction
p_shap <- function(data)
{
  tmp <- summary(data)
  tmp2 <- shapiro.test(tmp$residuals)
}

## résultats
df_shap <- df_stat |>
  split(df_stat$phylum) |> 
  map(fit_model) |>
  map(p_shap) |> 
  map_dfr(tidy, .id = "phylum") |> 
  dplyr::select(-method) |> 
  rename(W = statistic,
         W_pval = p.value)
df_shap
```

```
## # A tibble: 7 × 3
##   phylum               W   W_pval
##   <chr>            <dbl>    <dbl>
## 1 Actinobacteriota 0.953 0.0948  
## 2 Bacteroidota     0.978 0.625   
## 3 Desulfobacterota 0.874 0.000363
## 4 Firm_bact_ratio  0.891 0.00103 
## 5 Firmicutes       0.941 0.0382  
## 6 Fusobacteriota   0.958 0.139   
## 7 Proteobacteria   0.923 0.00942
```

**Lettres**


```r
## fonction
p_letters <- function(data)
{
  tmp <- emmeans(data, pairwise ~ Segment * Age)
  tmp2 <- cld(tmp$emmeans, alpha = 0.05, Letters = letters, adjust ="tukey")
}

## résultats
result_letters <- df_stat |>
  split(df_stat$phylum) |> 
  map(fit_model) |>
  map(p_letters) |> 
  map_dfr(tidy, .id = "phylum")
result_letters
```

```
## # A tibble: 28 × 9
##    phylum           Segment Age   estimate std.e…¹    df conf.low conf.…² .group
##    <chr>            <chr>   <chr>    <dbl>   <dbl> <dbl>    <dbl>   <dbl> <chr> 
##  1 Actinobacteriota Jejunum D35     0.192   0.0229  6.88  0.115    0.268  " a"  
##  2 Actinobacteriota Colon   D35     0.222   0.0354 15.0   0.121    0.322  " a"  
##  3 Actinobacteriota Jejunum D21     0.223   0.0229  6.74  0.146    0.300  " a"  
##  4 Actinobacteriota Colon   D21     0.243   0.0238  7.58  0.165    0.320  " a"  
##  5 Bacteroidota     Jejunum D35     0.175   0.0464  5.22  0.00242  0.347  " a " 
##  6 Bacteroidota     Jejunum D21     0.219   0.0464  5.17  0.0460   0.392  " a " 
##  7 Bacteroidota     Colon   D35     0.698   0.0562  9.37  0.526    0.871  "  b" 
##  8 Bacteroidota     Colon   D21     0.762   0.0470  5.44  0.590    0.933  "  b" 
##  9 Desulfobacterota Jejunum D21     0.0384  0.0200 10.6  -0.0216   0.0984 " a " 
## 10 Desulfobacterota Jejunum D35     0.0386  0.0200 10.7  -0.0213   0.0985 " a " 
## # … with 18 more rows, and abbreviated variable names ¹​std.error, ²​conf.high
```

**Mise en forme des résultats des lettres en colonnes**


```r
df_letters <- result_letters |>
  dplyr::select(phylum, Segment, Age, .group) |> 
  pivot_wider(names_from = c(Segment, Age), 
              values_from = .group)
df_letters
```

```
## # A tibble: 7 × 5
##   phylum           Jejunum_D35 Colon_D35 Jejunum_D21 Colon_D21
##   <chr>            <chr>       <chr>     <chr>       <chr>    
## 1 Actinobacteriota " a"        " a"      " a"        " a"     
## 2 Bacteroidota     " a "       "  b"     " a "       "  b"    
## 3 Desulfobacterota " a "       "  b"     " a "       "  b"    
## 4 Firm_bact_ratio  "  b"       " a "     "  b"       " a "    
## 5 Firmicutes       "  b d"     " a c "   "   cd"     " ab  "  
## 6 Fusobacteriota   " a "       " a "     " a "       "  b"    
## 7 Proteobacteria   " a"        " a"      " a"        " a"
```

### **Merge all results**


```r
left_join(df_shap, df_pval) |> 
  left_join(df_letters) -> full_result
full_result
```

```
## # A tibble: 7 × 14
##   phylum          W  W_pval (Inter…¹   Segment     Age Segme…² p_adj_S…³ p_adj…⁴
##   <chr>       <dbl>   <dbl>    <dbl>     <dbl>   <dbl>   <dbl>     <dbl>   <dbl>
## 1 Actinobact… 0.953 9.48e-2 8.68e-52 2.04e-  1 0.370    0.799  2.38e-  1  0.675 
## 2 Bacteroido… 0.978 6.25e-1 1.52e-50 2.23e-119 0.386    0.679  1.56e-118  0.675 
## 3 Desulfobac… 0.874 3.63e-4 2.48e-56 2.38e- 38 0.788    0.775  8.33e- 38  0.788 
## 4 Firm_bact_… 0.891 1.03e-3 1.31e-24 5.42e- 22 0.563    0.715  1.27e- 21  0.760 
## 5 Firmicutes  0.941 3.82e-2 0        2.22e- 14 0.254    0.264  3.89e- 14  0.675 
## 6 Fusobacter… 0.958 1.39e-1 1.24e- 6 9.45e-  2 0.00169  0.0302 1.32e-  1  0.0118
## 7 Proteobact… 0.923 9.42e-3 7.01e-18 9.81e-  1 0.652    0.979  9.81e-  1  0.760 
## # … with 5 more variables: p_adj_Age_Segment <dbl>, Jejunum_D35 <chr>,
## #   Colon_D35 <chr>, Jejunum_D21 <chr>, Colon_D21 <chr>, and abbreviated
## #   variable names ¹​`(Intercept)`, ²​`Segment:Age`, ³​p_adj_Segment, ⁴​p_adj_Age
```

```r
full_result |> 
  mutate(across(where(is.double), ~round(.x,4))) |> 
  datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-d9cd19af78ac36d9f2c3" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-d9cd19af78ac36d9f2c3">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7"],["Actinobacteriota","Bacteroidota","Desulfobacterota","Firm_bact_ratio","Firmicutes","Fusobacteriota","Proteobacteria"],[0.9528,0.9783,0.874,0.8906,0.9413,0.9577,0.9228],[0.0948,0.625,0.0004,0.001,0.0382,0.1392,0.0094],[0,0,0,0,0,0,0],[0.2037,0,0,0,0,0.0945,0.9812],[0.3703,0.3858,0.7882,0.5628,0.2539,0.0017,0.6517],[0.7995,0.6785,0.7748,0.7151,0.2635,0.0302,0.9788],[0.2376,0,0,0,0,0.1323,0.9812],[0.6752,0.6752,0.7882,0.7603,0.6752,0.0118,0.7603],[0.9327,0.9327,0.9327,0.9327,0.9223,0.2113,0.9788],[" a"," a "," a ","  b","  b d"," a "," a"],[" a","  b","  b"," a "," a c "," a "," a"],[" a"," a "," a ","  b","   cd"," a "," a"],[" a","  b","  b"," a "," ab  ","  b"," a"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>phylum<\/th>\n      <th>W<\/th>\n      <th>W_pval<\/th>\n      <th>(Intercept)<\/th>\n      <th>Segment<\/th>\n      <th>Age<\/th>\n      <th>Segment:Age<\/th>\n      <th>p_adj_Segment<\/th>\n      <th>p_adj_Age<\/th>\n      <th>p_adj_Age_Segment<\/th>\n      <th>Jejunum_D35<\/th>\n      <th>Colon_D35<\/th>\n      <th>Jejunum_D21<\/th>\n      <th>Colon_D21<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6,7,8,9,10]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```r
full_result |> 
  mutate(across(where(is.double), ~round(.x,4))) |> 
  flextable() -> tab 
  theme_vanilla(tab)
```

```{=html}
<template id="77bd2c9c-d2e9-4218-89c3-77f36797390d"><style>
.tabwid table{
  border-spacing:0px !important;
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  border-color: transparent;
  caption-side: top;
}
.tabwid-caption-bottom table{
  caption-side: bottom;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td, .tabwid th {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
.katex-display {
    margin: 0 0 !important;
}
</style><div class="tabwid"><style>.cl-f4be936a{}.cl-f4b6dd6e{font-family:'Arial';font-size:11pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-f4b6dd78{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-f4b9d000{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-f4b9d014{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-f4b9e482{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-f4b9e483{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-f4b9e48c{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-f4b9e496{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-f4b9e497{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-f4b9e4a0{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 0.5pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-f4b9e4a1{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-f4b9e4aa{width:0.75in;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0.5pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-f4be936a'><thead><tr style="overflow-wrap:break-word;"><th class="cl-f4b9e482"><p class="cl-f4b9d000"><span class="cl-f4b6dd6e">phylum</span></p></th><th class="cl-f4b9e483"><p class="cl-f4b9d014"><span class="cl-f4b6dd6e">W</span></p></th><th class="cl-f4b9e483"><p class="cl-f4b9d014"><span class="cl-f4b6dd6e">W_pval</span></p></th><th class="cl-f4b9e483"><p class="cl-f4b9d014"><span class="cl-f4b6dd6e">(Intercept)</span></p></th><th class="cl-f4b9e483"><p class="cl-f4b9d014"><span class="cl-f4b6dd6e">Segment</span></p></th><th class="cl-f4b9e483"><p class="cl-f4b9d014"><span class="cl-f4b6dd6e">Age</span></p></th><th class="cl-f4b9e483"><p class="cl-f4b9d014"><span class="cl-f4b6dd6e">Segment:Age</span></p></th><th class="cl-f4b9e483"><p class="cl-f4b9d014"><span class="cl-f4b6dd6e">p_adj_Segment</span></p></th><th class="cl-f4b9e483"><p class="cl-f4b9d014"><span class="cl-f4b6dd6e">p_adj_Age</span></p></th><th class="cl-f4b9e483"><p class="cl-f4b9d014"><span class="cl-f4b6dd6e">p_adj_Age_Segment</span></p></th><th class="cl-f4b9e482"><p class="cl-f4b9d000"><span class="cl-f4b6dd6e">Jejunum_D35</span></p></th><th class="cl-f4b9e482"><p class="cl-f4b9d000"><span class="cl-f4b6dd6e">Colon_D35</span></p></th><th class="cl-f4b9e482"><p class="cl-f4b9d000"><span class="cl-f4b6dd6e">Jejunum_D21</span></p></th><th class="cl-f4b9e482"><p class="cl-f4b9d000"><span class="cl-f4b6dd6e">Colon_D21</span></p></th></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-f4b9e48c"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">Actinobacteriota</span></p></td><td class="cl-f4b9e496"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9528</span></p></td><td class="cl-f4b9e496"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0948</span></p></td><td class="cl-f4b9e496"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0</span></p></td><td class="cl-f4b9e496"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.2037</span></p></td><td class="cl-f4b9e496"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.3703</span></p></td><td class="cl-f4b9e496"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.7995</span></p></td><td class="cl-f4b9e496"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.2376</span></p></td><td class="cl-f4b9e496"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.6752</span></p></td><td class="cl-f4b9e496"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9327</span></p></td><td class="cl-f4b9e48c"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a</span></p></td><td class="cl-f4b9e48c"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a</span></p></td><td class="cl-f4b9e48c"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a</span></p></td><td class="cl-f4b9e48c"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">Bacteroidota</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9783</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.6250</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0000</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.3858</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.6785</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0000</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.6752</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9327</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a </span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">  b</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a </span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">  b</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">Desulfobacterota</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.8740</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0004</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0000</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.7882</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.7748</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0000</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.7882</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9327</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a </span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">  b</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a </span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">  b</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">Firm_bact_ratio</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.8906</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0010</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0000</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.5628</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.7151</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0000</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.7603</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9327</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">  b</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a </span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">  b</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a </span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">Firmicutes</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9413</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0382</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0000</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.2539</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.2635</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0000</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.6752</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9223</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">  b d</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a c </span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">   cd</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> ab  </span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">Fusobacteriota</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9577</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.1392</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0945</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0017</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0302</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.1323</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0118</span></p></td><td class="cl-f4b9e4a0"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.2113</span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a </span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a </span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a </span></p></td><td class="cl-f4b9e497"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">  b</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-f4b9e4a1"><p class="cl-f4b9d000"><span class="cl-f4b6dd78">Proteobacteria</span></p></td><td class="cl-f4b9e4aa"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9228</span></p></td><td class="cl-f4b9e4aa"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.0094</span></p></td><td class="cl-f4b9e4aa"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0</span></p></td><td class="cl-f4b9e4aa"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9812</span></p></td><td class="cl-f4b9e4aa"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.6517</span></p></td><td class="cl-f4b9e4aa"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9788</span></p></td><td class="cl-f4b9e4aa"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9812</span></p></td><td class="cl-f4b9e4aa"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.7603</span></p></td><td class="cl-f4b9e4aa"><p class="cl-f4b9d014"><span class="cl-f4b6dd78">0.9788</span></p></td><td class="cl-f4b9e4a1"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a</span></p></td><td class="cl-f4b9e4a1"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a</span></p></td><td class="cl-f4b9e4a1"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a</span></p></td><td class="cl-f4b9e4a1"><p class="cl-f4b9d000"><span class="cl-f4b6dd78"> a</span></p></td></tr></tbody></table></div></template>
<div class="flextable-shadow-host" id="ea8757fb-65e3-47df-b9d8-17b6f6d64cfc"></div>
<script>
var dest = document.getElementById("ea8757fb-65e3-47df-b9d8-17b6f6d64cfc");
var template = document.getElementById("77bd2c9c-d2e9-4218-89c3-77f36797390d");
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```

**Save in excel sheet**


```r
write.xlsx(full_result,
           file = "result_table.xlsx",
           sheetName = "Stat model")
```

## **Remarque**

**Attention à l'ordre des lignes quand on fusionne des data.frame, éviter les `cbind()`**  

Utiliser la fonction `merge()` ou les fonctions `_join()` de dplyr qui garantissent le respect de l'ordre des lignes. La fusion s'effectue à partir d'une colonne commune aux 2 tables.

Exemple fusion des taxon et de la table d'abondance


```r
TA <- data.frame(otu_table(relatabunphy))
TAXA <- tax_table(tab_Phylum)[,"Phylum"]
```

**merge()**


```r
## merge()
temp_merge <- merge(TAXA , TA, by = "row.names")
```

**left_join()**

N'accepte que des data.frames
On utilise la fonction rownames_to_column() du package tibble qui transforme les rownames en colonne avec "rowname" en en-tête par défaut.
left_join peut repèrer automatiquement la colonne en commun: "rowname"  


```r
## left_join()
temp_join <- left_join(rownames_to_column(data.frame(TAXA)), rownames_to_column(TA))
```


```r
## par défaut merge() tri sur la colonne en commun
temp_merge[1:5,1:5]
```

```
##     Row.names           Phylum GPS158690_CCTTGA.JLGMN_L001
## 1   Cluster_1       Firmicutes                0.9984080886
## 2 Cluster_189 Campylobacterota                0.0000000000
## 3 Cluster_294  Patescibacteria                0.0000000000
## 4  Cluster_31     Bacteroidota                0.0002796601
## 5  Cluster_57   Fusobacteriota                0.0000000000
##   GPS158692_TCGTTC.JLGMN_L001 GPS158694.PCR450.6E_GTTTCT.JLGMN_L001
## 1                7.877842e-01                            0.88554166
## 2                9.093665e-05                            0.00000000
## 3                0.000000e+00                            0.00000000
## 4                1.965444e-01                            0.08483679
## 5                0.000000e+00                            0.00000000
```

```r
temp_join[1:5,1:5]
```

```
##       rowname           Phylum GPS158690_CCTTGA.JLGMN_L001
## 1  Cluster_31     Bacteroidota                0.0002796601
## 2 Cluster_189 Campylobacterota                0.0000000000
## 3   Cluster_1       Firmicutes                0.9984080886
## 4  Cluster_97 Desulfobacterota                0.0000000000
## 5  Cluster_86 Actinobacteriota                0.0004087340
##   GPS158692_TCGTTC.JLGMN_L001 GPS158694.PCR450.6E_GTTTCT.JLGMN_L001
## 1                1.965444e-01                            0.08483679
## 2                9.093665e-05                            0.00000000
## 3                7.877842e-01                            0.88554166
## 4                3.485905e-03                            0.00268524
## 5                8.487420e-04                            0.01476882
```

**Session Info**


```r
sessionInfo()
```

```
## R version 4.2.0 (2022-04-22 ucrt)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 19044)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=French_France.utf8  LC_CTYPE=French_France.utf8   
## [3] LC_MONETARY=French_France.utf8 LC_NUMERIC=C                  
## [5] LC_TIME=French_France.utf8    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] flextable_0.8.3     DT_0.26             emmeans_1.8.3      
##  [4] lmerTest_3.1-3      lme4_1.1-31         Matrix_1.5-3       
##  [7] multcomp_1.4-20     TH.data_1.1-1       MASS_7.3-58.1      
## [10] survival_3.4-0      mvtnorm_1.1-3       car_3.1-1          
## [13] carData_3.0-5       phyloseq_1.40.0     broom.mixed_0.2.9.4
## [16] broom_1.0.2         purrr_1.0.0         tidyr_1.2.1        
## [19] dplyr_1.0.10        tibble_3.1.8       
## 
## loaded via a namespace (and not attached):
##   [1] minqa_1.2.5            colorspace_2.0-3       ellipsis_0.3.2        
##   [4] estimability_1.4.1     XVector_0.36.0         base64enc_0.1-3       
##   [7] rstudioapi_0.14        listenv_0.9.0          furrr_0.3.1           
##  [10] fansi_1.0.3            xml2_1.3.3             codetools_0.2-18      
##  [13] splines_4.2.0          cachem_1.0.6           knitr_1.41            
##  [16] ade4_1.7-20            jsonlite_1.8.4         nloptr_2.0.3          
##  [19] pbkrtest_0.5.1         cluster_2.1.4          compiler_4.2.0        
##  [22] backports_1.4.1        assertthat_0.2.1       fastmap_1.1.0         
##  [25] cli_3.5.0              htmltools_0.5.4        tools_4.2.0           
##  [28] igraph_1.3.5           coda_0.19-4            gtable_0.3.1          
##  [31] glue_1.6.2             GenomeInfoDbData_1.2.8 reshape2_1.4.4        
##  [34] Rcpp_1.0.9             Biobase_2.56.0         jquerylib_0.1.4       
##  [37] vctrs_0.5.1            Biostrings_2.64.1      rhdf5filters_1.8.0    
##  [40] multtest_2.52.0        ape_5.6-2              nlme_3.1-161          
##  [43] crosstalk_1.2.0        iterators_1.0.14       xfun_0.36             
##  [46] stringr_1.5.0          globals_0.16.2         lifecycle_1.0.3       
##  [49] future_1.30.0          zlibbioc_1.42.0        zoo_1.8-11            
##  [52] scales_1.2.1           parallel_4.2.0         biomformat_1.24.0     
##  [55] sandwich_3.0-2         rhdf5_2.40.0           yaml_2.3.6            
##  [58] ggplot2_3.4.0          gdtools_0.2.4          sass_0.4.4            
##  [61] stringi_1.7.8          S4Vectors_0.34.0       foreach_1.5.2         
##  [64] permute_0.9-7          BiocGenerics_0.42.0    zip_2.2.2             
##  [67] boot_1.3-28.1          GenomeInfoDb_1.32.4    rlang_1.0.6           
##  [70] pkgconfig_2.0.3        systemfonts_1.0.4      bitops_1.0-7          
##  [73] evaluate_0.19          lattice_0.20-45        Rhdf5lib_1.18.2       
##  [76] htmlwidgets_1.6.0      tidyselect_1.2.0       parallelly_1.33.0     
##  [79] plyr_1.8.8             magrittr_2.0.3         R6_2.5.1              
##  [82] multcompView_0.1-8     IRanges_2.30.1         generics_0.1.3        
##  [85] DBI_1.1.3              withr_2.5.0            pillar_1.8.1          
##  [88] mgcv_1.8-41            abind_1.4-5            RCurl_1.98-1.9        
##  [91] crayon_1.5.2           uuid_1.1-0             utf8_1.2.2            
##  [94] rmarkdown_2.19         officer_0.5.0          grid_4.2.0            
##  [97] data.table_1.14.6      vegan_2.6-4            forcats_0.5.2         
## [100] digest_0.6.31          xtable_1.8-4           numDeriv_2016.8-1.1   
## [103] stats4_4.2.0           munsell_0.5.0          bslib_0.4.2
```


