# JetEngineAnalysis

## Purpose
This package was designed to meet and exceed the requirements of Georgia Tech's AE 4451 group project in Fall 2018. This package is capable of component by component thermodynamic analysis of a turbojet or turbofan engine. It can take any number of design inputs and produce the maximum specific thrust for that configuration. Alternatively, you can select three design requirements (specific thrust, altitude, and velocity) and the code will determine the desired engine component configuration and mass ratios (including compressor bleed). 


## Guide:

### Installation

* Download repository .zip to desired path on local machine.
* Unzip package to desired location
* Open Matlab to this location
* Navigate to the 'Project Files' folder

### Run the code
* Double click the 'GUI.mlapp' file in Matlab's 'Current Folder' display
  * This will open the graphic user interface
* In the app that opens, you will notice 4 tabs:
  #### Initialization:
    * These are the most common elements that will need to be changed
    * The values are self explanatory, but ensure you change everything necessary before executing
    * Based on selection of 'Engine type' and 'Bypass nozzle', several values will automatically populate
  #### Constants:
    * This tab contains values that may not necessarily need to be changed that often
    * Several of these may not be known by the designer beforehand and will be approximated
  #### Optimization:
    * Here you can visualize the effect varying the flight velocity has on the specific thrust and thrust specific fuel consumption
    * You can also use the knobs to change design parameters to see the effect this has on the specific thrust.
    * Here you can also select a desired specific thrust, velocity, and altitude; the code will automaticall find the optimal mass ratios and display the results on the graph. 
  #### Results:
    * This tab shows all the calculated intermediary and final results in a single table.
### Closing the app
* When you run the code, an excel file named "Results.xlsx" will be written to your local machine in the project's file path.
* When you close the app you will be prompted if you wish to delete this file.



