

export FuzCav=../Packages/FuzCav 
echo "FuzCav Dir:"
echo $FuzCav 
count=0
for i in $1/*/*site.mol2  ; do
    echo "$i"
    sitefile="$i"
    $FuzCav/utils/CaTagger.pl $sitefile > tmp/site1_Tagged.mol2 
    echo "tmp/site1_Tagged.mol2"  > "tmp/listCavTagged"
    dir=`dirname "$i"`
    java -jar $FuzCav/dist/3pointPharCav.jar -d $FuzCav/utils/resDef/tableDefCA.txt -t $FuzCav/utils/triplCav/interval.txt -l tmp/listCavTagged -o $dir/FPCount.txt -c
    cat tmp/listCavTagged 
    rm tmp/listCavTagged
done 

