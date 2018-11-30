# JetEngineAnalysis

## Purpose
This package was designed to meet and exceed the requirements of Georgia Tech's AE 4451 group project in the Fall 2018 semester. The code is capable of component by component thermodynamic analysis of a turbojet, turbofan, or ramjet engine using **MatLab version 2018b**. 

Out code can take any number of design inputs and determine the maximum specific thrust and minimum thrust specific fuel consumption for that configuration. Alternatively, you can select your desired specific thrust, altitude, and mach number and the code will calculate the desired engine component configuration and mass ratios (including compressor bleed) that will minimize the thrust specific fuel consumption. 

*This application requires the aerospace toolbox be installed in MatLab. Please see [this](https://www.mathworks.com/help/aerotbx/) link for information on the toolbox and [this](https://www.mathworks.com/matlabcentral/answers/101885-how-do-i-install-additional-toolboxes-into-an-existing-installation-of-matlab) link for how to install it.*
___
## Tutorial

[![Youtube Tutorial](https://img.youtube.com/vi/c-5Q0Uuz6wk/maxresdefault.jpg)](http://www.youtube.com/watch?v=c-5Q0Uuz6wk)
___
## How it Works:

  Every file you need to run the application is located in the "Project Files" folder within the repository. Everything is controlled by the GUI.mlapp file. This application will take user inputs and call the necessary functions. The JetPro_Project function will take perform the bulk of the calculations. The outputs of each engine component (implemented as local functions in the code) become the inputs for the next component.
  
The "optimization" function takes user inputs, varies six components, and finds their values that minimize Thrust Specific Fuel Consumption (TSFC) while meeting two user-defined flight conditions. It uses MatLab's GlobalSearch function to call fmincon (using a sequential quadratic programming 'sqp' algorithm) with design related constraints and determines what engine configuration (afterburner or no) meets these constraints. 
    
Similarly, the "Maximize" function will keep physical parameters constant (such as bypass ratio and fan pressure ratio) and vary three components to maximize the Specific Thrust (ST) still subject to design constraints. This is again accomplished via GlobalSearch calling fmincon which in turn calls the JetPro_Project function.
    
*Note: GlobalSearch uses a scatter search algorithm to generate a number of start points, allowing us to obtain a global minimum.*
___
## Guide:

### Installation

* Download repository .zip to desired path on local machine.
* Unzip resulting package to desired location.
* Open Matlab to this location.
* Navigate to the 'Project Files' folder.

### Run the code
* Double click the 'GUI.mlapp' file in Matlab's 'Current Folder'. After several seconds this will open the graphic user interface.
* In the app that opens, you will notice 4 tabs:
  #### Initialization:
    * These are the most common elements that will need to be changed.
    * The values are self explanatory, but ensure you change everything necessary before executing.
    * Based on selected aircraft from dropdown menu, several values will automatically populate.
    * By clicking 'Calculate' the code will execute and populate the results throughout the app.
    * By clicking 'Optimize' the code will find the minimum TSFC that meets two user-defined flight conditions. Results will populate the respective fields, the visualization tab, and the table in the 'Results' tab.
    * By clicking the 'Maximize' button, combustor fuel-air ratio, afterburner fuel-air ratio, and bleed ratio are varied to maximize the specific thrust. This will also update the visualization and results tabs.
  #### Constants:
    * This tab contains values that may not necessarily need to be changed that often.
    * Several of these may not be known by the designer beforehand and the default values can be used.
  #### Visualization:
    * Here you can visualize the effect varying the fuel-air ratio, afterburner fuel-air ratio, and bleed ratio has on the specific thrust and thrust specific fuel consumption at different Mach numbers.
    * The two figures on the right show the temperature and pressure at the exit of each component.
    * The  'Optimize' button on the 'Initialization' tab will also update these graphs.
  #### Results:
    * This tab shows all the calculated intermediary and final results in a single table. Essentially the raw data.
### Closing the app
* When you run the code, an excel file named "Results.xlsx" will be written to your local machine in MatLab's current working directory. When you close the app you will be prompted if you wish to delete this file.



