library(officer)

doc <- read_docx()
print(doc,target="resources/templates/r_template.docx")
