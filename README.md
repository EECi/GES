# GES
Greenhouse Energy Simulation

## Introduction
This MATLAB code simulates heat, mass and CO<sub>2</sub> exchange in an unheated, ventilated single-zone greenhouse for a simple test case.  The heat transfer processes simulated include convection, conduction and radiation, together with plant transpiration for simulation of heat and mass exchange due to evaporation of water from the leaf surface.  Simple models of photosynthesis and crop growth for tomatoes are included in order to simulate CO<sub>2</sub> exchange.  The model is based on the GDGCM (Pieters, J. and Deltour, J., 'The Gembloux Dynamic Greenhouse Climate Model - GDGCM', 2000 ) and on the thesis by Vanthoor (Vanthoor, B.H.,'A model-based greenhouse design method', PhD Thesis, Wageningen University, 2011).

## How to run model
The files included here are as follows:

* GESModel.m   

  This is the main MATLAB file which initialises and calls the differential equation solver ode15s with derivatives.m.  Using parameters.m and the weather data in SampleWeather.csv, running GESModel.m provides time histories of temperature, air moisture content and CO<sub>2</sub> concentration. 

* derivatives.m

  Contains the governing heat, mass and CO<sub>2</sub> balance equations for solution at each timestep
  
* parameters.m

  Contains parameters for greenhouse and plant geometry and material properties 

* climterp_linear.m

  Sub-routine called by GESModel.m to interpolate weather data according to the timestep specified in parameters.m

* lamorturb.m

  Sub-routine called by derivatives.m in order to calculate whether flow is laminar or turbulent for convection calculation

* sat_conc.m

  Sub-routine which converts relative humidity into air moisture content.

* PlotResults.m

  Routine for plotting time histories of temperature, moisture content, relative humidity and CO<sub>2</sub> concentration output by GESModel.m

* SampleWeather.csv

  Hourly input weather data, in the format 'Hour, Ambient Temperature (<sup>o</sup>C), Sky Temperature (<sup>o</sup>C), Windspeed (m/s), Relative Humidity (%), Direct Solar Radiation (order NE wall, NE Roof, SE Wall, SE Roof, SW Wall, SW Roof, NW Wall, NW Roof) (W/m<sup>2</sup>), Diffuse Solar Radiation (order as for Direct)(W/m<sup>2</sup>)
  
The model has been run with MATLAB R2018a https://uk.mathworks.com/products/matlab.html
 
## Sample data
The example given is a single zone unheated greenhouse, dimensions 10 x 25 x 5m with no artificial lighting ior supplemental CO<sub>2</sub>, as illustrated below.![alt Greenhouse](https://github.com/EECi/GES/blob/master/images/Greenhouse.PNG)

The heat, mass and CO<sub>2</sub> exchanges are modelled as indicated below.![alt Balance](https://github.com/EECi/GES/blob/master/images/Balance.PNG)


 ## Project status
 
 This code is under development and has been modified for use in multi-zone, rooftop and underground greenhouses.  More information can be found under the University of Cambridge Engineering Department Energy Efficient Cities initative website here: https://www.eeci.cam.ac.uk/research/greenhouse-energy-simulation
 
