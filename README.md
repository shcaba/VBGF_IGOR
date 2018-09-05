# IGOR+: Fitting Growth Curves with Random Effects

IGOR is a shiny tool that allows you to do the following tasks:
1. Create an Age vs. Otolith Weight model (piecewise or linear) with user data,
2. Use an existing Age vs. Otolith Weight model to predict ages for fish with known otolith weights,
3. Run a growth curve model (gompertz, linear, logistic, schnute, or von bertlanffy) to fit fish lengths with ages.

## Prerequisite

This Shiny app requires users to have a C++ compiler installed.

## How To Run The App



## How Do I Run an Age vs. Otolith Weight Analysis?
<p float = "left">
  <img src = "imgs/1.png" width = "420" height = "400" />
  <img src = "imgs/2.png" width = "420" height = "400" /> 
</p>
When you start this app, you should see two panels for file upload. To run an Age vs. Otolith Weight model, you should upload a file with specified format in the panel shown above. Then in the Age vs. Otolith Weight Data tab, you would be able to see a control panel on the left hand side that provides you options for the model and your input data table on the right hand side.
<br></br>
After you choose the starting values for the model, you would see the result of the run in the Summary & Plot tab. You can also download the model for later use. In addition, the predicted ages will appear in the Fitted Values tab and you could this data for a growth curve analysis.
<p float = "left">
  <img src = "imgs/3.png" width = "420" height = "400" />
  <img src = "imgs/4.png" width = "420" height = "400" /> 
</p>

## How Do I Use An Existing Age vs. Otolith Weight Model To Predict Ages?

You now have an Age vs. Otolith Model and a data file that has otolith weights and you want to predict the ages for these fish. In the Use An Existing Otolith Weight vs. Age Model tab, you should see a file upload panel where you would upload an R object of type `lm` and a csv file that contains all your data. Your predicted ages will appear in a data table on the right.
<p align = "center">
  <img src = "imgs/5.png" />
</p>

## How Do I Run A Growth Curve Model Analysis?

There are two ways to feed in fish age and length data: 
1. There are two file upload panels when you start the app and to run a growth curve model analysis. You can upload a file with specified format as Length vs. Age data.
2. If you have run an Age vs. Otolith Weight model, you can use the predicted data as fish age and length data for this analysis.
<br></br>
The raw data from upload or predictions from otolith weight will show up in the Input Data tab under the Length vs. Age Data tab. You can then filter the data as you like and confirm the data choice. The filtered data will show up in the Selected Data tab under the Length vs. Age Data tab.
<p float="left">
  <img src = "imgs/6.png" />
  <img src = "imgs/7.png" /> 
</p>
After you have selected your data, you can run a growth curve analysis with a model you like. A scatterplot and a fitting curve will show up on the right. The screenshot below shows a standard run with Von Bertlanffy Model. You can also run a model with random effects.
<p align = "center">
  <img src = "imgs/8.png" />
</p>
You can find out the results (parameter estimates) of your model runs in the Growth Curve Summaries tab.
<p align = "center">
  <img src = "imgs/9.png" />
</p>
