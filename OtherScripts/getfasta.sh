echo -e "24296" | while read G; do curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&id=${G}" | grep -A 1 "<Link>" | grep "<Id>" | cut -d '>' -f 2 | cut -d '<' -f 1 | while read S ; do curl -s "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${S}&retmode=text&rettype=fasta" ; done;  done