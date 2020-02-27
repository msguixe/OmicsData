from Bio import Entrez

Entrez.email = "your.email@here.com"

with Entrez.esearch(db="pubmed", term="gwas") as search_handle:
    search_data = Entrez.read(search_handle)
list_of_ids = search_data["IdList"]

with Entrez.efetch(db="pubmed", id=list_of_ids, retmode="xml") as fetch_handle:
    data = Entrez.read(fetch_handle)

articles = data["PubmedArticle"]
authors, years, titles = [], [], []
for art in articles:
    authors.append(art["MedlineCitation"]["Article"]["AuthorList"][0]["LastName"])
    titles.append(art["MedlineCitation"]["Article"]["ArticleTitle"])
    if len(art["MedlineCitation"]["Article"]["ArticleDate"]) > 0:
        years.append(art["MedlineCitation"]["Article"]["ArticleDate"][0]["Year"])
    else:
        years.append("XXXX")

for i, title in enumerate(titles):
    # print(i, authors[i], years[i], title)
    print(f"{authors[i]:15} {years[i]}   {title}")
