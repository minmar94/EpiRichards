# EpiRichards
Modelling epidemic outbreaks using the Richards' curve. In particular, Monkeypox, COVID-19, Campylobacteriosis and Measles. 

## Data
The $\texttt{data}$ folder includes:

    - covid19 positive cases for Italian Omicron BA.5 wave
    - adjacency matrix of Italian regions base on geographical proximity
    - Italian regional population on the log-scale scaled by a factor of $10^{-4}$
    - measles infections for Germany from 2005 to 2018 at the federal state level
    - monkeypox cases for Italy and United States from mid May to mid September 2022
    
## Fit
The $\texttt{data}$ folder includes 4 subfolders, each one for each application to a different epidemic with the codes for fitting a GLM model with the Richards' curve
