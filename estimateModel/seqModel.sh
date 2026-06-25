Folder="xx/"
echo $Folder
if [ ! -f "$Folder/finalModEsts.rds" ]; then
echo "Starting model"
echo "Stage 1"
#RunRIID.R
nice Rscript RunR.R $Folder
fi
echo "Stage 1 complete"
if [ ! -f "$Folder/modEsts.rds" ]; then
nice Rscript convertData.R $Folder
fi
echo "Data made"
echo "Stage 2"
nice Rscript RF.R $Folder
echo "Stage 2 complete"