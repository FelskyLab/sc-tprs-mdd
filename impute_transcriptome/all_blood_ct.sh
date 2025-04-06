count=0
for filename in ~/sc-tprs-mdd/immune_cell_models/*.db; do
  fn=`basename $filename .db`
  echo $fn
  count=$((count+1))
done
echo $count