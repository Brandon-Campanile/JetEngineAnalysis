# JetEngineAnalysis

## Purpose
This package was designed to meet and exceed the requirements of Georgia Tech's AE 4451 group project in Fall 2018. This package is capable of component by component thermodynamic analysis of a turbojet or turbofan engine. It can take any number of design inputs and determine the maximum specific thrust and minimum thrust specific fuel consumption for that configuration. Alternatively, you can select your desired specific thrust, altitude, and mach number and the code will calculate the desired engine component configuration and mass ratios (including compressor bleed) that will minimize the thrust specific fuel consumption. 


## Guide:

### Installation

* Download repository .zip to desired path on local machine.
* Unzip resulting package to desired location.
* Open Matlab to this location.
* Navigate to the 'Project Files' folder.

### Run the code
* Double click the 'GUI.mlapp' file in Matlab's 'Current Folder'. This will open the graphic user interface.

* In the app that opens, you will notice 4 tabs:
  #### Initialization:
    * These are the most common elements that will need to be changed.
    * The values are self explanatory, but ensure you change everything necessary before executing.
    * Based on selection of 'Engine type' and 'Bypass nozzle', several values will automatically populate.
    * By clicking 'Calculate' the code will run and the table in the 'Results' tab will populate.
    * By clicking 'Optimize' the code will find the minimum TSFC that meets the input ST. Results will populate the respective field and table in the 'Results' tab. *Note this will take several minutes.
  #### Constants:
    * This tab contains values that may not necessarily need to be changed that often.
    * Several of these may not be known by the designer beforehand and the default values can be used.
  #### Visualization:
    * Here you can visualize the effect varying the fuel-air ratio, afterburner fuel-air ratio, and bleed ratio has on the specific thrust and thrust specific fuel consumption.
    * The  'Optimize' button on the 'Initialization' tab will also update these graphs.
  #### Results:
    * This tab shows all the calculated intermediary and final results in a single table.
### Closing the app
* When you run the code, an excel file named "Results.xlsx" will be written to your local machine in the project's file path.
* When you close the app you will be prompted if you wish to delete this file.



