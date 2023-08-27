using Documenter, combo_seq, DocumenterTools

makedocs(sitename="NEXTFLEX Combo-Seq Documentation")

deploydocs(
    repo = "github.com/cookienocreams/combo_seq.git"
    , branch = "docs"
    , devbranch = "main"
)