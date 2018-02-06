for FILE in `ls ../ModellerRSC758/Data/*pdb`;do 
    echo $FILE
    freesasa --radii naccess --depth atom $FILE --format pdb -o Data/`basename $FILE`.sasa;
done
