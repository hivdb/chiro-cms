#! /bin/sh

cd $(dirname $0)/..
DRDB_VERSION=$(\grep drdbVersion pages/sierra-sars2.yml | awk '{print $2}' | jq -r | sed 's/-slim$//g')
mkdir -p downloads/resistance-mutations

rm -rf local/covid-drdb.db
curl --fail -sSLo local/covid-drdb.db "https://s3-us-west-2.amazonaws.com/cms.hivdb.org/covid-drdb/covid-drdb-${DRDB_VERSION}.db" --compressed
rm -rf downloads/resistance-mutations/latest.tsv
cat <<EOF | sqlite3 local/covid-drdb.db
.headers on
.mode tabs
.output downloads/resistance-mutations/latest.tsv
SELECT
  im.gene,
  im.position,
  CASE im.amino_acid
    WHEN 'ins' THEN '_'
    WHEN 'del' THEN '-'
    WHEN 'stop' THEN '*'
    ELSE im.amino_acid
  END AS aa
FROM isolate_mutations im
  JOIN isolate_pairs ip ON
    im.iso_name = ip.iso_name
  JOIN isolates iso ON
    im.iso_name = iso.iso_name
  JOIN susc_results s ON
    s.control_iso_name=ip.control_iso_name AND
    s.iso_name = ip.iso_name
WHERE
  (
    (
      ip.num_mutations = 1 AND
      fold >= 5
    ) OR
    EXISTS (
      SELECT 1
      FROM antibody_epitopes abe, antibodies ab
      WHERE
        abe.ab_name=ab.ab_name AND
        ab.visibility=1 AND
        abe.position=im.position
    )
  ) AND
  EXISTS (
    SELECT 1
    FROM rx_antibodies rxab, antibodies ab
    WHERE
      s.ref_name=rxab.ref_name AND
      s.rx_name=rxab.rx_name AND
      rxab.ab_name=ab.ab_name AND
      visibility = 1
  ) AND
  im.gene = 'S' AND
  (
    iso.var_name IS NULL OR 
    iso.var_name NOT IN (
      'SARS-CoV-1', 'WIV1', 'pCoV-GD', 'pCoV-GX', 'bCoV-RaTG13', 'Unknown Variant'
    )
  ) AND
  NOT EXISTS (
    SELECT 1
    FROM ignore_mutations igm
    WHERE
      igm.gene=im.gene AND
      igm.position=im.position AND
      igm.amino_acid=im.amino_acid
  )
GROUP BY im.gene, im.position, aa
ORDER BY im.gene, im.position, aa;
.quit
EOF
echo "Create: $(realpath downloads/resistance-mutations)/latest.tsv"

cat <<EOF | python3
import csv
import json

with open('downloads/resistance-mutations/latest.tsv') as fp:
    rows = []
    for row in csv.DictReader(fp, delimiter='\t'):
        row['position'] = int(row['position'])
        rows.append(row)

with open('downloads/resistance-mutations/latest.json', 'w') as fp:
    json.dump({'MAB': rows}, fp, indent=2)

EOF
echo "Create: $(realpath downloads/resistance-mutations)/latest.json"
