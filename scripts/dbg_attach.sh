EXEC=./partition_crs

bat process.txt

proc_id=`tail -1 process.txt`
lldb -p $proc_id $EXEC 