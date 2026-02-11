Projekt na poprawę oceny z zajęć „Proteomika” w zakresie sprawdzenia umiejętności przeliczania 
zestawu danych proteomicznych. 
Zadaniem Studenta jest opracowanie danych z eksperymentu proteomicznego zdeponowanego w 
repozytorium PRIDE o nr. referencyjnym PXD000279. Jest to zestaw danych użytych w pracy Cox et 
al. z 2014 roku (https://doi.org/10.1074/mcp.M113.031591). Interesuje nas jedynie część danych 
użytych w eksperymencie „Dynamic Range Benchmark”. 
Cel:  
Z danych MaxQuant (proteinGroups) wykonać analizę różnicową UPS2 vs UPS1, a potem policzyć 
empiryczny FDR na białkach E. coli i dokładność fold-change dla UPS. 
UPS1 to rekombinowane białka ludzkie w mieszaninie ekwimolarnej, a UPS2 to białka w grupach o 
znanych proporcjach w zakresie ~5 rzędów wielkości 
(https://www.sigmaaldrich.com/PL/pl/technical-documents/technical-article/protein
biology/protein-mass-spectrometry/ups1-and-ups2
proteomic?srsltid=AfmBOoomK_wlyK5JN2ynCk_mnrcB1TQpnqkIVIis-EBIUmrgWm8_qaG7) 
Tło E. coli w porównaniu UPS2 vs UPS1 powinno być ~stałe → a więc stanowić idealne “false positives 
control”. 
Plan pracy: 
1) Dane: pobierz z FTP projektu PXD000279 pliki MaxQuant output (proteinGroups.txt + ew. 
experimentalDesignTemplate.txt). Link do FTP jest na stronie datasetu. Interesuje Cię część 
UPS1/UPS2 spiked in E. coli (dynamic range benchmark). 
2) Preprocessing w R lub Perseus: wczytanie proteinGroups.txt, filtrowanie, log2 transform, 
filtrowanie danych numerycznych “min. 2 z 3 powtórzeń” (albo podobny), ew. normalizacja (np. 
median lub vsn), imputacja albo analiza bez imputacji (wariant porównawczy) 
3) Różnicowanie UPS2 vs UPS1: statystyka, wykresy: volcano, heatmap dla top białek, PCA po 
normalizacji  
4) Benchmark: empiryczny FDR na tle E. coli (oznacz białka jako E. coli vs UPS (Homo sapiens) (często 
w identyfikatorach MaxQuant/UniProt widać _ECOLI / _HUMAN; jeśli nie, to mapowanie po 
FASTA/ID), policz: ile białek E. coli wyszło jako istotne przy FDR<0.05 to daje nam empiryczny odsetek 
fałszywych trafień, dokładność fold-change dla UPS, dla białek UPS porównaj obserwowany log2FC z 
oczekiwanym (UPS2 ma znane grupy stężeń względem UPS1, informacje na podanej wcześniej 
stronie producenta). 
metryki: korelacja, RMSE, wykres “expected vs observed”. 
5) [AI / SOTA] Analiza szumu tła E. coli: Użycie modelu językowego białek (ESM-2) do wygenerowania 
embeddingów dla białek E. coli. Wizualizacja przestrzeni latentnej (UMAP/PCA) w celu sprawdzenia, 
czy "fałszywie pozytywne" białka tła grupują się biologicznie (np. czy są to białka rybosomalne, 
błonowe, czy o specyficznej strukturze). To pokazuje nowoczesne podejście "AI-driven biology".
Dodatkowo może być „benchmark” decyzji pipeline’u (najlepiej 2–3 warianty), np. porównaj: median 
normalization + bez imputacji, vsn + QRILC (lub MinProb), wskaż który wariant daje mniej False 
Positive w E. coli i lepszą zgodność FC dla UPS. 
Raport: 
Powinien zawierać opis danych, opis „benchmarku”, rysunki (pełna dowolność, np. „missingness 
heatmap”, PCA, volcano (kolory: E. coli vs UPS), expected-vs-observed, tabela false positive), krótka 
sekcja “Wnioski”: który wariant pipeline’u jest najbardziej wiarygodny i dlaczego. 
Status:
[X] 1. Pobranie danych
[X] 2. Preprocessing
[X] 3. Analiza różnicowa
[X] 4. Benchmark
[X] 5. Analiza AI (AI Latent Space)