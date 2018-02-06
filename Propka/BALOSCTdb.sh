for FILE in `ls ../ModellerBALOSCTdb/Data/*pdb`;do 
    echo $FILE
    /root/bin/propka31 -q $FILE 
    mv `basename $FILE|sed -e 's/pdb/pka/'` Data/
    rm *propka_input
done
