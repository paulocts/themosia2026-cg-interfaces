xmax=$(awk '!/^[@#]/{if($1>m)m=$1} END{print m}' density.xvg)
# 2a) header (all @ and # lines)
awk '/^[@#]/{print}' density.xvg > density_unwrapped.xvg

# 2b) process + unwrap data only, then sort numerically by x
awk -v xmax="$xmax" '
!/^[@#]/{
  x=$1
  if (x < -2.0) x = x + 2*xmax
  printf "%10.5f", x
  for(i=2;i<=NF;i++) printf " %12.6f", $i
  printf "\n"
}' density.xvg | sort -n >> density_unwrapped.xvg
mv density_unwrapped.xvg density.xvg
