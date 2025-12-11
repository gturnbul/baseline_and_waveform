1. Ensure that the makefile has the correct program.

2. $ make
   Using the make file, this will create a file called OffsetCalc.exe which can then be run. 

3. $ ./ OffsetCalc.exe #run
  This will create offset CSV files for Memory cell, Mod16 and end spike calculations, with labels;

  run_#runnumber_EoW_offsets.csv, 
  run_#runnumber_Kept_events_per_om_bin.csv, 
  run_#runnumber_MemCell_offsets.csv, 
  run_#runnumber_Mod16_offsets.csv

To apply MemCell offsets, you need to reorder the waveform to memory cell order using the FCR before applying the offsets. If using mod16 offsets, time order is sufficient.
