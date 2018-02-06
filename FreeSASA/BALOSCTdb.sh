for FILE in `ls ../ModellerBALOSCTdb/Data/*pdb`;do 
    echo $FILE
    #freesasa --radii naccess --depth atom $FILE --format pdb -o Data/`basename $FILE`.sasa;
    freesasa --depth atom $FILE --format pdb -o Data/`basename $FILE`.sasa;
done
