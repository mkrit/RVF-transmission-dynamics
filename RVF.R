library(shiny)
require(deSolve)
#--------------- Call Cpp functions
Rcpp::sourceCpp('CppFunctions.cpp')
require(Rcpp)
options(warn=-1)

# Define UI for application that draws a histogram
ui = navbarPage("Rift Valley Fever",
       tabPanel("Run simulation & select plots",
         mainPanel(fluidRow(br()),
           fluidRow(column(12, actionButton("runIt", label = h2("Run simulation")),
                    br(), br(),
                    sliderInput("year", label = "Number of years to run simulation",
                                min = 1, max = 50, value = 27, step = 1), align="center")),
           fluidRow(column(4, numericInput("plotStart", h4("Start plotting at (day):"), value = 6841)), 
                    column(4, numericInput("plotEnd", h4("Stop plotting at (day):") , value = 9720)),
                    column(4, numericInput("rngSeed", h4("RNG seed"), value = 123))),
           fluidRow(br(), hr()),
           fluidRow(column(3, h1("species"), hr(), h1("human"), hr(), h1("animal"),hr(),
                           h1("vector A"), hr(), h1("vector B"), hr(), h1("vector C"), hr(), h1("vector D", hr()),
                           align="center"),
                    column(3, h1("zone 1"), hr(),
                           h1(checkboxInput("plotH1", NULL, value = TRUE)),  hr(),
                           h1(checkboxInput("plotM1", NULL, value = TRUE)),  hr(),
                           h1(checkboxInput("plotA1", NULL, value = TRUE)),  hr(),
                           h1(checkboxInput("plotB1", NULL, value = FALSE)), hr(),
                           h1(checkboxInput("plotC1", NULL, value = TRUE)),  hr(),
                           h1(checkboxInput("plotD1", NULL, value = FALSE)), hr(),
                           align="center"),
                    column(3, h1("zone 2"), hr(),
                           h1(checkboxInput("plotH2", NULL, value = TRUE)),  hr(),
                           h1(checkboxInput("plotM2", NULL, value = TRUE)),  hr(),
                           h1(checkboxInput("plotA2", NULL, value = FALSE)), hr(),
                           h1(checkboxInput("plotB2", NULL, value = TRUE)),  hr(),
                           h1(checkboxInput("plotC2", NULL, value = FALSE)), hr(),
                           h1(checkboxInput("plotD2", NULL, value = TRUE)),  hr(),
                           align="center"),
                    column(3, h1("zone 3"), hr(),
                           h1(checkboxInput("plotH3", NULL, value = TRUE)),  hr(),
                           h1(checkboxInput("plotM3", NULL, value = TRUE)),  hr(),
                           h1(checkboxInput("plotA3", NULL, value = FALSE)), hr(),
                           h1(checkboxInput("plotB3", NULL, value = TRUE)),  hr(),
                           h1(checkboxInput("plotC3", NULL, value = FALSE)), hr(),
                           h1(checkboxInput("plotD3", NULL, value = TRUE)),
                           align="center")
             ) # fluidRow
           ) # mainPanel
         ), # tabPanel plot
       tabPanel("Initial state",
         mainPanel(
           fluidRow(
                    column(3, h3("Compartment"),
                           hr(), h3("HS", style="padding:4px;"), hr(), h3("HE", style="padding:4px;"),
                           hr(), h3("HI", style="padding:4px;"), hr(), h3("HR", style="padding:4px;"), 
                           hr(), h3("MS", style="padding:4px;"), hr(), h3("ME", style="padding:4px;"),
                           hr(), h3("MI", style="padding:4px;"), hr(), h3("MR", style="padding:4px;"),
                           hr(), h3("AQ", style="padding:4px;"), hr(), h3("AP", style="padding:4px;"),
                           hr(), h3("AS", style="padding:4px;"), hr(), h3("AI", style="padding:4px;"),
                           hr(), h3("BQ", style="padding:4px;"), hr(), h3("BP", style="padding:4px;"),
                           hr(), h3("BS", style="padding:4px;"), hr(), h3("BI", style="padding:4px;"),
                           hr(), h3("CP", style="padding:4px;"), hr(), h3("CS", style="padding:4px;"),
                           hr(), h3("CI", style="padding:4px;"),
                           hr(), h3("DP", style="padding:4px;"), hr(), h3("DS", style="padding:4px;"),
                           hr(), h3("DI", style="padding:4px;"), hr(), align="center"), # column labels
                    column(3, h3("Zone 1"), hr(),
                           numericInput("HS1", NULL, value = 0), hr(), numericInput("HE1", NULL, value = 0), hr(),
                           numericInput("HI1", NULL, value = 0), hr(), numericInput("HR1", NULL, value = 0), hr(),
                           numericInput("MS1", NULL, value = 0), hr(), numericInput("ME1", NULL, value = 0), hr(),
                           numericInput("MI1", NULL, value = 0), hr(), numericInput("MR1", NULL, value = 0), hr(),
                           numericInput("AQ1", NULL, value = 100), hr(), numericInput("AP1", NULL, value = 9900), hr(),
                           numericInput("AS1", NULL, value = 0), hr(), numericInput("AI1", NULL, value = 0), hr(),
                           numericInput("BQ1", NULL, value = 0), hr(), numericInput("BP1", NULL, value = 0), hr(),
                           numericInput("BS1", NULL, value = 0), hr(), numericInput("BI1", NULL, value = 0), hr(),
                           numericInput("CP1", NULL, value = 100), hr(), numericInput("CS1", NULL, value = 0), hr(),
                           numericInput("CI1", NULL, value = 0), hr(),
                           numericInput("DP1", NULL, value = 0), hr(), numericInput("DS1", NULL, value = 0), hr(),
                           numericInput("DI1", NULL, value = 0), hr(),
                           align="center"), # column zone 1
                    column(3, h3("Zone 2"), hr(),
                           numericInput("HS2", NULL, value = 1000), hr(), numericInput("HE2", NULL, value = 0), hr(),
                           numericInput("HI2", NULL, value = 0), hr(), numericInput("HR2", NULL, value = 0), hr(),
                           numericInput("MS2", NULL, value = 2500), hr(), numericInput("ME2", NULL, value = 0), hr(),
                           numericInput("MI2", NULL, value = 0), hr(), numericInput("MR2", NULL, value = 0), hr(),
                           numericInput("AQ2", NULL, value = 0), hr(), numericInput("AP2", NULL, value = 0), hr(),
                           numericInput("AS2", NULL, value = 0), hr(), numericInput("AI2", NULL, value = 0), hr(),
                           numericInput("BQ2", NULL, value = 0), hr(), numericInput("BP2", NULL, value = 10), hr(),
                           numericInput("BS2", NULL, value = 0), hr(), numericInput("BI2", NULL, value = 0), hr(),
                           numericInput("CP2", NULL, value = 0), hr(), numericInput("CS2", NULL, value = 0), hr(),
                           numericInput("CI2", NULL, value = 0), hr(),
                           numericInput("DP2", NULL, value = 0), hr(), numericInput("DS2", NULL, value = 1000), hr(),
                           numericInput("DI2", NULL, value = 0), hr(),
                           align="center"), # column zone 2
                    column(3, h3("Zone 3"), hr(),
                           numericInput("HS3", NULL, value = 0), hr(), numericInput("HE3", NULL, value = 0), hr(),
                           numericInput("HI3", NULL, value = 0), hr(), numericInput("HR3", NULL, value = 0), hr(),
                           numericInput("MS3", NULL, value = 0), hr(), numericInput("ME3", NULL, value = 0), hr(),
                           numericInput("MI3", NULL, value = 0), hr(), numericInput("MR3", NULL, value = 0), hr(),
                           numericInput("AQ3", NULL, value = 0), hr(), numericInput("AP3", NULL, value = 0), hr(),
                           numericInput("AS3", NULL, value = 0), hr(), numericInput("AI3", NULL, value = 0), hr(),
                           numericInput("BQ3", NULL, value = 0), hr(), numericInput("BP3", NULL, value = 10), hr(),
                           numericInput("BS3", NULL, value = 0), hr(), numericInput("BI3", NULL, value = 0), hr(),
                           numericInput("CP3", NULL, value = 0), hr(), numericInput("CS3", NULL, value = 0), hr(),
                           numericInput("CI3", NULL, value = 0), hr(),
                           numericInput("DP3", NULL, value = 1000), hr(), numericInput("DS3", NULL, value = 0), hr(),
                           numericInput("DI3", NULL, value = 0), hr(),
                           align="center") # column zone 3
             ) # fluidRow
           ) # mainPanel
         ), # tabPanel state
       tabPanel("General model information",
         mainPanel(
           fluidRow(column(12, checkboxInput("elNino", h4("El Nino flooding"), value = TRUE, width="800px"), align="left")),
           fluidRow(
                    column(4, sliderInput("d5", h4("Start"), min = 1, max = 360, value = 270), align="left"),
                    column(4, sliderInput("d6", h4("End"), min = 1, max = 360, value = 315), align="left")),
           fluidRow(column(12, checkboxInput("flood", h4("Annual flooding"), value = TRUE, width="800px"), align="left")),
           fluidRow(column(4, sliderInput("d3", h4("Start"), min = 1, max = 360, value = 91), align="left"),
                    column(4, sliderInput("d4", h4("End"), min = 1, max = 360, value = 120), align="left"),
                    column(4, numericInput("flood_prop", h4("Proportion"), value = 0.05/30))),
           fluidRow(column(12, checkboxInput("seasonHatch", h4("Seasonal effect hatching"), value = TRUE, width="800px"), align="left")),
           fluidRow(column(4, sliderInput("ds", h4("Delay"), min = -180, max = 180, value = 0), align="left"),
                    column(4, sliderInput("nPeak", h4("Number of annual peaks"), min = 1, max = 4, value = 1), align="left")),
           fluidRow(column(12, checkboxInput("wetDry", h4("Annual variation climate"), value = TRUE, width="800px"), align="left")),
           fluidRow(column(4, sliderInput("w", h4("Dry years"), min = 0, max = 5, value = 3), align="left"),
                    column(4, sliderInput("W", h4("Total years"), min = 0, max = 10, value = 7), align="left"),
                    column(4, sliderInput("mmm", h4("Minimum"), min = 0, max = 1, value = 0.1))),
           fluidRow(column(12, checkboxInput("transHumance", h4("Annual transhumance"), value = TRUE, width="800px"), align="left")),
           fluidRow(column(4, sliderInput("d1", h4("Homestead to pasture"), min = 1, max = 360, value = 181), align="left"),
                    column(4, sliderInput("d2", h4("Pasture to homestead"), min = 1, max = 360, value = 330), align="left")),
           fluidRow(column(12, checkboxInput("shearing", h4("Increased susceptibility of animals"), value = FALSE, width="800px"), align="left")),
           fluidRow(column(4, sliderInput("shearBeg", h4("Start"), min = 1, max = 360, value = 1), align="left"),
                    column(4, sliderInput("shearEnd", h4("End"), min = 1, max = 360, value = 1), align="left"),
                    column(4, numericInput("shearUp", h4("Factor"), value = 1))),
           fluidRow(column(12, numericInput("b_wl", h4("Infection rate wildlife"), value = 0))),
           fluidRow(column(12, numericInput("O_alt", h4("Number of bites on alternative hosts"), value = 0)))
           ) # mainPanel
         ), # tabPanel general
       tabPanel("People",
         mainPanel(
           fluidRow(column(6, numericInput("g_h", "Birth rate", value = 0.0001),
                           numericInput("m_h", "Natural mortality rate", value = 0.0001),
                           hr(),
                           numericInput("x_h", "Average length incubation period", value = 4),
                           numericInput("a_h", "Average length infective period", value = 3),
                           numericInput("d_h", "Disease-specific mortality rate", value = 0.01),
                           numericInput("r_h", "Average length immune period", value = 900),
                           hr(),
                           numericInput("l_h12", "Migration rate zones 1->2", value = 0.05),
                           numericInput("l_h13", "Migration rate zones 1->3", value = 0.0001),
                           numericInput("l_h21", "Migration rate zones 2->1", value = 0.05),
                           numericInput("l_h23", "Migration rate zones 2->3", value = 0.001),
                           numericInput("l_h31", "Migration rate zones 3->1", value = 0.005),
                           numericInput("l_h32", "Migration rate zones 3->2", value = 0.05),
                           hr(),
                           align="left"), # first column
                    column(6, numericInput("p_ha", "Infection transfer rate human -> vector A", value = 0.89),
                           numericInput("p_hb", "Infection transfer rate human -> vector B", value = 0.89),
                           numericInput("p_hc", "Infection transfer rate human -> vector C", value = 0.81),
                           numericInput("p_hd", "Infection transfer rate human -> vector D", value = 0.81),
                           hr(),
                           numericInput("f_mh1", "Contact rate human-animal zone 1", value = 2.5),
                           numericInput("f_mh2", "Contact rate human-animal zone 2", value = 2.5),
                           numericInput("f_mh3", "Contact rate human-animal zone 3", value = 2.5),
                           hr(),
                           numericInput("h_h1", "Maximum supported biting rate zone 1", value = 25),
                           numericInput("h_h2", "Maximum supported biting rate zone 2", value = 25),
                           numericInput("h_h3", "Maximum supported biting rate zone 3", value = 25),
                           hr(),
                           align="left") # second column
             ) # fluidRow
           ) # mainPanel
         ), # tabPanel human
       tabPanel("Animal",
         mainPanel(
           fluidRow(column(6, numericInput("g_m_u", "Birth rate non-infected", value = 0.00082),
                           numericInput("p_a_i", "Proportion abortion due to RVF", value = 0.9),
                           numericInput("m_m", "Natural mortality rate", value = 0.0008),
                           hr(),
                           numericInput("k_m1", "Carrying capacity zone 1", value = 500000),
                           numericInput("k_m2", "Carrying capacity zone 2", value = 500000),
                           numericInput("k_m3", "Carrying capacity zone 3", value = 500000),
                           hr(),
                           numericInput("x_m", "Average length incubation period", value = 24/3.25),
                           numericInput("a_m", "Average length infective period", value = 20),
                           numericInput("d_m", "Disease-specific mortality rate", value = 0.05),
                           numericInput("r_m", "Average length immune period", value = 900),
                           hr(),
                           numericInput("l_m12Base", "Migration rate zones 1->2*", value = 1e-5),
                           numericInput("l_m21Base", "Migration rate zones 2->1*", value = 1e-5),
                           numericInput("l_m13", "Migration rate zones 1->3", value = 0),
                           numericInput("l_m23", "Migration rate zones 2->3", value = 0.0001),
                           numericInput("l_m31", "Migration rate zones 3->1", value = 0),
                           numericInput("l_m32", "Migration rate zones 3->2", value = 0.005), hr(),
                           align="left"), # first column,
                    column(6, numericInput("p_mh00", "Infection transfer rate animal -> human", value = 0.001),
                           numericInput("p_ma", "Infection transfer rate animal -> vector A", value = 0.89),
                           numericInput("p_mb", "Infection transfer rate animal -> vector B", value = 0.89),
                           numericInput("p_mc", "Infection transfer rate animal -> vector C", value = 0.81),
                           numericInput("p_md", "Infection transfer rate animal -> vector D", value = 0.81),
                           hr(),
                           numericInput("h_m", "Maximum supported biting rate", value = 50),
                           hr(),
                           align="left")
             ) # fluidRow
           ) # mainPanel
         ), # tabPanel animal
       tabPanel("Vector A",
         mainPanel(
           fluidRow(column(6, numericInput("g_a", "Oviposition rate", value = 10),
                    numericInput("z_a", "Proportion vertical transmission", value = 0.5),
                    numericInput("t_a", "Hatching rate", value = 0.2),
                    hr(),
                    numericInput("m_a", "Mortality rate adults", value = 1/3),
                    numericInput("m_aq1", "Mortality rate infected eggs zone 1", value = 0.00001),
                    numericInput("m_aq2", "Mortality rate infected eggs zone 2", value = 0.00001),
                    numericInput("m_aq3", "Mortality rate infected eggs zone 3", value = 0.00001),
                    numericInput("m_ap1", "Mortality rate uninfected eggs zone 1", value = 0.00001),
                    numericInput("m_ap2", "Mortality rate uninfected eggs zone 2", value = 0.00001),
                    numericInput("m_ap3", "Mortality rate uninfected eggs zone 3", value = 0.00001),
                    hr(),
                    numericInput("v_a", "Maximum biting rate", value = 0.5),
                    numericInput("e_ah", "Proportion feeding on humans", value = 0.1),
                    numericInput("p_ah", "Infection transfer rate vector A -> human", value = 0.01),
                    numericInput("e_am", "Proportion feeding on animals", value = 0.3),
                    numericInput("p_am", "Infection transfer rate vector A -> animal", value = 0.01),
                    hr(),
                    align="left"), # first column
             column(6, numericInput("k_a1", "Carrying capacity zone 1", value = 175000),
                    numericInput("k_a2", "Carrying capacity zone 2", value = 175000),
                    numericInput("k_a3", "Carrying capacity zone 3", value = 175000),
                    hr(),
                    numericInput("l_a12", "Migration rate zones 1->2", value = 0),
                    numericInput("l_a13", "Migration rate zones 1->3", value = 0),
                    numericInput("l_a21", "Migration rate zones 2->1", value = 0),
                    numericInput("l_a23", "Migration rate zones 2->3", value = 0),
                    numericInput("l_a31", "Migration rate zones 3->1", value = 0),
                    numericInput("l_a32", "Migration rate zones 3->2", value = 0),
                    hr(),
                    align="left") # second column
             ) # fluidRow
           ) # mainPanel
         ), # tabPanel Vector A
       tabPanel("Vector B",
         mainPanel(
           fluidRow(column(6, numericInput("g_b", "Oviposition rate", value = 25),
                    numericInput("z_b", "Proportion vertical transmission", value = 0.05),
                    numericInput("t_b", "Hatching rate", value = 0.2),
                    hr(),
                    numericInput("m_b", "Mortality rate adults", value = 0.1),
                    numericInput("m_bq1", "Mortality rate infected eggs zone 1", value = 0.005),
                    numericInput("m_bq2", "Mortality rate infected eggs zone 2", value = 0.005),
                    numericInput("m_bq3", "Mortality rate infected eggs zone 3", value = 0.005),
                    numericInput("m_bp1", "Mortality rate uninfected eggs zone 1", value = 0.005),
                    numericInput("m_bp2", "Mortality rate uninfected eggs zone 2", value = 0.005),
                    numericInput("m_bp3", "Mortality rate uninfected eggs zone 3", value = 0.005),
                    hr(),
                    numericInput("v_b", "Maximum biting rate", value = 0.5),
                    numericInput("e_bh", "Proportion feeding on humans", value = 0.01),
                    numericInput("p_bh", "Infection transfer rate vector B -> human", value = 0.01),
                    numericInput("e_bm", "Proportion feeding on animals", value = 0.25),
                    numericInput("p_bm", "Infection transfer rate vector B -> animal", value = 0.01),
                    hr(),
                    align="left"), # first column
             column(6, numericInput("k_b1", "Carrying capacity zone 1", value = 175000),
                    numericInput("k_b2", "Carrying capacity zone 2", value = 175000),
                    numericInput("k_b3", "Carrying capacity zone 3", value = 175000),
                    hr(),
                    numericInput("l_b12", "Migration rate zones 1->2", value = 0),
                    numericInput("l_b13", "Migration rate zones 1->3", value = 0),
                    numericInput("l_b21", "Migration rate zones 2->1", value = 0),
                    numericInput("l_b23", "Migration rate zones 2->3", value = 0),
                    numericInput("l_b31", "Migration rate zones 3->1", value = 0),
                    numericInput("l_b32", "Migration rate zones 3->2", value = 0),
                    hr(),
                    align="left") # second column
             ) # fluidRow
           ) # mainPanel
         ), # tabPanel Vector B
       tabPanel("Vector C",
         mainPanel(
           fluidRow(column(6,
                    numericInput("g_c", "Oviposition rate", value = 25),
                    numericInput("t_c", "Hatching rate", value = 0.2),
                    hr(),
                    numericInput("m_c", "Mortality rate adults", value = 0.1),
                    numericInput("m_cp1", "Mortality rate eggs zone 1", value = 0.002),
                    numericInput("m_cp2", "Mortality rate eggs zone 2", value = 0.002),
                    numericInput("m_cp3", "Mortality rate eggs zone 3", value = 0.002),
                    hr(),
                    numericInput("v_c", "Maximum biting rate", value = 1),
                    numericInput("e_ch", "Proportion feeding on humans", value = 0.0025),
                    numericInput("p_ch", "Infection transfer rate vector C -> human", value = 0.07),
                    numericInput("e_cm", "Proportion feeding on animals", value = 0.02),
                    numericInput("p_cm", "Infection transfer rate vector C -> animal", value = 0.07),
                    hr(),
                    align="left"), # first column
             column(6, numericInput("k_c1", "Carrying capacity zone 1", value = 1750),
                    numericInput("k_c2", "Carrying capacity zone 2", value = 1750),
                    numericInput("k_c3", "Carrying capacity zone 3", value = 1750),
                    hr(),
                    numericInput("l_c12", "Migration rate zones 1->2", value = 0),
                    numericInput("l_c13", "Migration rate zones 1->3", value = 0),
                    numericInput("l_c21", "Migration rate zones 2->1", value = 0),
                    numericInput("l_c23", "Migration rate zones 2->3", value = 0),
                    numericInput("l_c31", "Migration rate zones 3->1", value = 0),
                    numericInput("l_c32", "Migration rate zones 3->2", value = 0),
                    hr(),
                    align="left") # second column
             ) # fluidRow
           ) # mainPanel
         ), # tabPanel Vector C
       tabPanel("Vector D",
         mainPanel(
           fluidRow(column(6, numericInput("g_d", "Oviposition rate", value = 25),
                    numericInput("t_d", "Hatching rate", value = 0.2),
                    hr(),
                    numericInput("m_d", "Mortality rate adults", value = 0.1),
                    numericInput("m_dp1", "Mortality rate eggs zone 1", value = 0.002),
                    numericInput("m_dp2", "Mortality rate eggs zone 2", value = 0.002),
                    numericInput("m_dp3", "Mortality rate eggs zone 3", value = 0.002),
                    hr(),
                    numericInput("v_d", "Maximum biting rate", value = 1),
                    numericInput("e_dh", "Proportion feeding on humans", value = 0.005),
                    numericInput("p_dh", "Infection transfer rate vector D -> human", value = 0.07),
                    numericInput("e_dm", "Proportion feeding on animals", value = 0.12),
                    numericInput("p_dm", "Infection transfer rate vector D -> animal", value = 0.07),
                    hr(),
                    align="left"), # first column
             column(6, numericInput("k_d1", "Carrying capacity zone 1", value = 17500),
                    numericInput("k_d2", "Carrying capacity zone 2", value = 17500),
                    numericInput("k_d3", "Carrying capacity zone 3", value = 17500),
                    hr(),
                    numericInput("l_d12", "Migration rate zones 1->2", value = 0),
                    numericInput("l_d13", "Migration rate zones 1->3", value = 0),
                    numericInput("l_d21", "Migration rate zones 2->1", value = 0),
                    numericInput("l_d23", "Migration rate zones 2->3", value = 0),
                    numericInput("l_d31", "Migration rate zones 3->1", value = 0),
                    numericInput("l_d32", "Migration rate zones 3->2", value = 0),
                    hr(),
                    align="left") # second column
             ) # fluidRow
           ) # mainPanel
         ), # tabPanel Vector D
       tabPanel("Summary tables",
         mainPanel(
           fluidRow(column(12, wellPanel(h4("Hosts"), tableOutput("summary"))), # column
                    column(12, wellPanel(h4("Vectors"), tableOutput("summary2"))), # column
                    column(12, wellPanel(h4("Seroprevalence"), tableOutput("summary3"))) # column
             ) # fluidRow
           ) # mainPanel
         ) # tabPanel Summary
  ) # navbarPage


# Define server logic required to draw a histogram
server = function(input, output) {
observe({ if(input$runIt > 0)
{
  start = Sys.time(); print(start)
  isolate({
    # Which plots to be generated
    ps = input$plotStart; pe = input$plotEnd; if(input$rngSeed>0) set.seed(input$rngSeed)
    plotH1 = input$plotH1; plotM1 = input$plotM1; plotA1 = input$plotA1; plotB1 = input$plotB1; plotC1 = input$plotC1; plotD1 = input$plotD1
    plotH2 = input$plotH2; plotM2 = input$plotM2; plotA2 = input$plotA2; plotB2 = input$plotB2; plotC2 = input$plotC2; plotD2 = input$plotD2
    plotH3 = input$plotH3; plotM3 = input$plotM3; plotA3 = input$plotA3; plotB3 = input$plotB3; plotC3 = input$plotC3; plotD3 = input$plotD3
    
    # Initialization of all the parameters and the initial values of the compartments per species
    aa = as.double(input$HS1); ab = as.double(input$HE1); ac = as.double(input$HI1); ad = as.double(input$HR1)
    ae = as.double(input$HS2); af = as.double(input$HE2); ag = as.double(input$HI2); ah = as.double(input$HR2)
    ai = as.double(input$HS3); aj = as.double(input$HE3); ak = as.double(input$HI3); al = as.double(input$HR3)
    am = as.double(input$MS1); an = as.double(input$ME1); ao = as.double(input$MI1); ap = as.double(input$MR1)
    aq = as.double(input$MS2); ar = as.double(input$ME2); as = as.double(input$MI2); at = as.double(input$MR2)
    au = as.double(input$MS3); av = as.double(input$ME3); aw = as.double(input$MI3); ax = as.double(input$MR3)
    ay = as.double(input$AQ1); az = as.double(input$AP1); ba = as.double(input$AS1); bb = as.double(input$AI1)
    bc = as.double(input$AQ2); bd = as.double(input$AP2); be = as.double(input$AS2); bf = as.double(input$AI2)
    bg = as.double(input$AQ3); bh = as.double(input$AP3); bi = as.double(input$AS3); bj = as.double(input$AI3)
    bk = as.double(input$BQ1); bl = as.double(input$BP1); bm = as.double(input$BS1); bn = as.double(input$BI1)
    bo = as.double(input$BQ2); bp = as.double(input$BP2); bq = as.double(input$BS2); br = as.double(input$BI2)
    bs = as.double(input$BQ3); bt = as.double(input$BP3); bu = as.double(input$BS3); bv = as.double(input$BI3)
    bw = as.double(input$CP1); bx = as.double(input$CS1); by = as.double(input$CI1)
    bz = as.double(input$CP2); ca = as.double(input$CS2); cb = as.double(input$CI2)
    cc = as.double(input$CP3); cd = as.double(input$CS3); ce = as.double(input$CI3)
    cf = as.double(input$DP1); cg = as.double(input$DS1); ch = as.double(input$DI1)
    ci = as.double(input$DP2); cj = as.double(input$DS2); ck = as.double(input$DI2)
    cl = as.double(input$DP3); cm = as.double(input$DS3); cn = as.double(input$DI3)
    # set ON/OFF switches and general parameters
    elNino = input$elNino; d5 = input$d5; d6 = input$d6
    flood = input$flood; d3 = input$d3; d4 = input$d4; flood_prop = input$flood_prop
    seasonHatch = input$seasonHatch; ds = input$ds; nPeak = 180/input$nPeak
    wetDry = input$wetDry; w = input$w;  W = input$W; mmm = input$mmm
    transHumance = input$transHumance; d1 = input$d1; d2 = input$d2
    shearing = input$shearing; shearBeg = input$shearBeg; shearEnd = input$shearEnd; shearUp = input$shearUp
    b_wl = input$b_wl; O_alt = input$O_alt; year = input$year
    if(ps < 1 | ps > year*360) ps = 1
    if(pe < ps | pe > year*360) pe = year*360
    ps = 1 + ps/0.1; pe = 1 + pe/0.1
    
    # Population parameter initialisations
    #---------------- Parameter initialisation human equations
    g_h = input$g_h; m_h = input$m_h; x_h = 1/input$x_h; a_h = (1-input$d_h)/input$a_h; d_h = input$d_h/input$a_h; p_mh00 = input$p_mh00
    f_mh1 = input$f_mh1; f_mh2 = input$f_mh2; f_mh3 = input$f_mh3; h_h1 = input$h_h1; h_h2 = input$h_h2; h_h3 = input$h_h3
    p_ha = input$p_ha; p_hb = input$p_hb; p_hc = input$p_hc; p_hd = input$p_hd; r_h = 1/input$r_h
    l_h12 = input$l_h12; l_h13 = input$l_h13; l_h21 = input$l_h21; l_h23 = input$l_h23; l_h31 = input$l_h31; l_h32 = input$l_h32
    #---------------- Parameter initialisation animal host equations
    g_m_u = input$g_m_u; g_m_i = (1-input$p_a_i)*g_m_u; m_m = input$m_m; x_m = input$x_m; a_m = (1-input$d_m)/input$a_m
    d_m = input$d_m/input$a_m; h_m = input$h_m
    p_ma = input$p_ma; p_mb = input$p_mb; p_mc = input$p_mc; p_md = input$p_md
    r_m = 1/input$r_m; k_m1 = input$k_m1; k_m2 = input$k_m2; k_m3 = input$k_m3
    l_m13 = input$l_m13; l_m23 = input$l_m23; l_m31 = input$l_m31; l_m32 = input$l_m32
    l_m12Base = input$l_m12Base; l_m21Base = input$l_m21Base
    #---------------- Parameter initialisation vector A equations
    g_a = input$g_a; z_a = input$z_a; m_a = input$m_a; t_a = input$t_a; v_a = input$v_a; e_ah = input$e_ah
    e_am = input$e_am; p_ah = input$p_ah; p_am = input$p_am
    k_a1 = input$k_a1; m_aq1 = input$m_aq1; m_ap1 = input$m_ap1; k_a2 = input$k_a2; m_aq2 = input$m_aq2; m_ap2 = input$m_ap2
    k_a3 = input$k_a3; m_aq3 = input$m_aq3; m_ap3 = input$m_ap3
    l_a12 = input$l_a12; l_a13 = input$l_a13; l_a21 = input$l_a21; l_a23 = input$l_a23; l_a31 = input$l_a31; l_a32 = input$l_a32
    #---------------- Parameter initialisation vector B equations
    g_b = input$g_b; z_b = input$z_b; m_b = input$m_b; t_b = input$t_b;  v_b = input$v_b
    e_bh = input$e_bh; e_bm = input$e_bm; p_bh = input$p_bh; p_bm = input$p_bm
    k_b1 = input$k_b1; m_bq1 = input$m_bq1; m_bp1 = input$m_bp1; k_b2 = input$k_b2; m_bq2 = input$m_bq2; m_bp2 = input$m_bp2
    k_b3 = input$k_b3; m_bq3 = input$m_bq3; m_bp3 = input$m_bp3
    l_b12 = input$l_b12; l_b13 = input$l_b13; l_b21 = input$l_b21; l_b23 = input$l_b23; l_b31 = input$l_b31; l_b32 = input$l_b32
    #---------------- Parameter initialisation vector C equations
    g_c = input$g_c; m_c = input$m_c; t_c = input$t_c;  v_c = input$v_c
    e_ch = input$e_ch; e_cm = input$e_cm; p_ch = input$p_ch; p_cm = input$p_cm
    k_c1 = input$k_c1; m_cp1 = input$m_cp1; k_c2 = input$k_c2; m_cp2 = input$m_cp2; k_c3 = input$k_c3; m_cp3 = input$m_cp3
    l_c12 = input$l_c12; l_c13 = input$l_c13; l_c21 = input$l_c21; l_c23 = input$l_c23; l_c31 = input$l_c31; l_c32 = input$l_c32
    #---------------- Parameter initialisation vector D equations
    g_d = input$g_d; m_d = input$m_d; t_d = input$t_d;  v_d = input$v_d
    e_dh = input$e_dh; e_dm = input$e_dm; p_dh = input$p_dh; p_dm = input$p_dm
    k_d1 = input$k_d1; m_dp1 = input$m_dp1; k_d2 = input$k_d2; m_dp2 = input$m_dp2; k_d3 = input$k_d3; m_dp3 = input$m_dp3
    l_d12 = input$l_d12; l_d13 = input$l_d13; l_d21 = input$l_d21; l_d23 = input$l_d23; l_d31 = input$l_d31; l_d32 = input$l_d32
  })
  
  # General initialisations
  #---------------- Initialisation of maximum rate allowed in ODE to avoid negative values
  #---------------- (if problems are encountered lower max_rate to 9)
  max_rate = 10
  #---------------- Timeframe
  times = round(seq(0, year*360, 0.1), 1); ntimes = length(times)
  #---------------- Wet and dry years
  dry2 = runif(year+1);  c0 = w/W
  #---------------- Initial state
  state = c(HS1 = aa, HE1 = ab, HI1 = ac, HR1 = ad, HS2 = ae, HE2 = af, HI2 = ag, HR2 = ah, HS3 = ai, HE3 = aj, HI3 = ak, HR3 = al,
            MS1 = am, ME1 = an, MI1 = ao, MR1 = ap, MS2 = aq, ME2 = ar, MI2 = as, MR2 = at, MS3 = au, ME3 = av, MI3 = aw, MR3 = ax,
            AQ1 = ay, AP1 = az, AS1 = ba, AI1 = bb, AQ2 = bc, AP2 = bd, AS2 = be, AI2 = bf, AQ3 = bg, AP3 = bh, AS3 = bi, AI3 = bj,
            BQ1 = bk, BP1 = bl, BS1 = bm, BI1 = bn, BQ2 = bo, BP2 = bp, BS2 = bq, BI2 = br, BQ3 = bs, BP3 = bt, BS3 = bu, BI3 = bv,
            CP1 = bw, CS1 = bx, CI1 = by,           CP2 = bz, CS2 = ca, CI2 = cb,           CP3 = cc, CS3 = cd, CI3 = ce,
            DP1 = cf, DS1 = cg, DI1 = ch,           DP2 = ci, DS2 = cj, DI2 = ck,           DP3 = cl, DS3 = cm, DI3 = cn)
  #---------------- Parameters to be passed to ODE function
  parameters = c(max_rate, flood_prop, b_wl, O_alt, d1, d2, d3, d4, d5, d6, year, c0, mmm, ds, nPeak, seasonHatch,
                 g_h, m_h, x_h, a_h, d_h, p_mh00, f_mh1, f_mh2, f_mh3, h_h1, h_h2, h_h3,
                 p_ha, p_hb, p_hc, p_hd, r_h, l_h12, l_h13, l_h21, l_h23, l_h31, l_h32,
                 g_m_u, g_m_i, m_m, x_m, a_m, d_m, h_m, p_ma, p_mb, p_mc, p_md, r_m, k_m1, k_m2, k_m3,
                 l_m13, l_m23, l_m31, l_m32,
                 g_a, z_a, m_a, v_a, e_ah, e_am, p_ah, p_am, k_a1, k_a2, k_a3, 
                 m_aq1, m_aq2, m_aq3, m_ap1, m_ap2, m_ap3, l_a12, l_a13, l_a21, l_a23, l_a31, l_a32,
                 g_b, z_b, m_b, t_b, v_b, e_bh, e_bm, p_bh, p_bm, k_b1, k_b2, k_b3,
                 m_bq1, m_bq2, m_bq3, m_bp1, m_bp2, m_bp3, l_b12, l_b13, l_b21, l_b23, l_b31, l_b32,
                 g_c, m_c, t_c, v_c, e_ch, e_cm, p_ch, p_cm, k_c1, k_c2, k_c3, m_cp1, m_cp2, m_cp3,
                 l_c12, l_c13, l_c21, l_c23, l_c31, l_c32,
                 g_d, m_d, t_d, v_d, e_dh, e_dm, p_dh, p_dm, k_d1, k_d2, k_d3, m_dp1, m_dp2, m_dp3,
                 l_d12, l_d13, l_d21, l_d23, l_d31, l_d32, shearing, shearBeg, shearEnd, shearUp,
                 wetDry, flood, elNino, transHumance, t_a,l_m21Base,l_m12Base, dry2)
  
  #--------------- ODE function
  rvfODE = function(t,  state, parameters)
  {
    with(as.list(c(state,parameters)),
         {ODE(t,  state, param=parameters)})
  }
  
  out = as.data.frame(ode(y=state, times=times, func=rvfODE, parms=parameters, method="rk4"))
  print(Sys.time() - start)
  
  # plot people in zone 1
  if(plotH1) drawPlot4(times[ps:pe], out$HS1[ps:pe], out$HR1[ps:pe], out$HE1[ps:pe], out$HI1[ps:pe],
                       "Human zone 1: susceptible recovered", "Exposed infected", c("Susceptible","Recovered","Exposed","Infected"), ps, pe, year)
  if(plotH2) drawPlot4(times[ps:pe], out$HS2[ps:pe], out$HR2[ps:pe], out$HE2[ps:pe], out$HI2[ps:pe],
                       "Human zone 2: susceptible recovered", "Exposed infected", c("Susceptible","Recovered","Exposed","Infected"), ps, pe, year)
  if(plotH3) drawPlot4(times[ps:pe], out$HS3[ps:pe], out$HR3[ps:pe], out$HE3[ps:pe], out$HI3[ps:pe],
                       "Human zone 3: susceptible recovered", "Exposed infected", c("Susceptible","Recovered","Exposed","Infected"), ps, pe, year)
  if(plotM1) drawPlot4(times[ps:pe], out$MS1[ps:pe], out$MR1[ps:pe], out$ME1[ps:pe], out$MI1[ps:pe],
                       "Animals zone 1: susceptible recovered", "Exposed infected", c("Susceptible","Recovered","Exposed","Infected"), ps, pe, year)
  if(plotM2) drawPlot4(times[ps:pe], out$MS2[ps:pe], out$MR2[ps:pe], out$ME2[ps:pe], out$MI2[ps:pe],
                       "Animals zone 2: susceptible recovered", "Exposed infected", c("Susceptible","Recovered","Exposed","Infected"), ps, pe, year)
  if(plotM2) drawPlot4(times[ps:pe], out$MS3[ps:pe], out$MR3[ps:pe], out$ME3[ps:pe], out$MI3[ps:pe],
                       "Animals zone 3: susceptible recovered", "Exposed infected", c("Susceptible","Recovered","Exposed","Infected"), ps, pe, year)
  if(plotA1) drawPlot4(times[ps:pe], out$AP1[ps:pe], out$AS1[ps:pe], out$AQ1[ps:pe], out$AI1[ps:pe],
                       "Vector A zone 1: uninfected eggs susceptible adults", "Infected eggs, infected adults",
                       c("Uninfected egg","Susceptible adult","Infected egg","Infected adult"), ps, pe, year)
  if(plotA2) drawPlot4(times[ps:pe], out$AP2[ps:pe], out$AS2[ps:pe], out$AQ2[ps:pe], out$AI2[ps:pe],
                       "Vector A zone 2: uninfected eggs susceptible adults", "Infected eggs, infected adults",
                       c("Uninfected egg","Susceptible adult","Infected egg","Infected adult"), ps, pe, year)
  if(plotA3) drawPlot4(times[ps:pe], out$AP3[ps:pe], out$AS3[ps:pe], out$AQ3[ps:pe], out$AI3[ps:pe],
                       "Vector A zone 3: uninfected eggs susceptible adults", "Infected eggs, infected adults",
                       c("Uninfected egg","Susceptible adult","Infected egg","Infected adult"), ps, pe, year)
  if(plotB1) drawPlot4(times[ps:pe], out$BP1[ps:pe], out$BS1[ps:pe], out$BQ1[ps:pe], out$BI1[ps:pe],
                       "Vector B zone 1: uninfected eggs susceptible adults", "Infected eggs, infected adults",
                       c("Uninfected egg","Susceptible adult","Infected egg","Infected adult"), ps, pe, year)
  if(plotB2) drawPlot4(times[ps:pe], out$BP2[ps:pe], out$BS2[ps:pe], out$BQ2[ps:pe], out$BI2[ps:pe],
                       "Vector B zone 2: uninfected eggs susceptible adults", "Infected eggs, infected adults",
                       c("Uninfected egg","Susceptible adult","Infected egg","Infected adult"), ps, pe, year)
  if(plotB3) drawPlot4(times[ps:pe], out$BP3[ps:pe], out$BS3[ps:pe], out$BQ3[ps:pe], out$BI3[ps:pe],
                       "Vector B zone 3: uninfected eggs susceptible adults", "Infected eggs, infected adults",
                       c("Uninfected egg","Susceptible adult","Infected egg","Infected adult"), ps, pe, year)
  if(plotC1) drawPlot3(times[ps:pe], out$CP1[ps:pe], out$CS1[ps:pe], out$CI1[ps:pe],
                       "Vector C zone 1: uninfected eggs and adults", "Infected adults", c("Egg", "Susceptible adult","Infected adult"), ps, pe, year)
  if(plotC2) drawPlot3(times[ps:pe], out$CP2[ps:pe], out$CS2[ps:pe], out$CI2[ps:pe],
                       "Vector C zone 2: uninfected eggs and adults", "Infected adults", c("Egg", "Susceptible adult","Infected adult"), ps, pe, year)
  if(plotC3) drawPlot3(times[ps:pe], out$CP3[ps:pe], out$CS3[ps:pe], out$CI3[ps:pe],
                       "Vector C zone 3: uninfected eggs and adults", "Infected adults", c("Egg", "Susceptible adult","Infected adult"), ps, pe, year)
  if(plotD1) drawPlot3(times[ps:pe], out$DP1[ps:pe], out$DS1[ps:pe], out$DI1[ps:pe],
                       "Vector D zone 1: uninfected eggs and adults", "Infected adults", c("Egg", "Susceptible adult","Infected adult"), ps, pe, year)
  if(plotD2) drawPlot3(times[ps:pe], out$DP2[ps:pe], out$DS2[ps:pe], out$DI2[ps:pe],
                       "Vector D zone 2: uninfected eggs and adults", "Infected adults", c("Egg", "Susceptible adult","Infected adult"), ps, pe, year)
  if(plotD3) drawPlot3(times[ps:pe], out$DP3[ps:pe], out$DS3[ps:pe], out$DI3[ps:pe],
                       "Vector D zone 3: uninfected eggs and adults", "Infected adults", c("Egg", "Susceptible adult","Infected adult"), ps, pe, year)
  ## Generate summarising tables
  ### Hosts
  summ = as.data.frame(array(NA, c(4,4), dimnames = list(c("Human", "", "Mammal", " "), c("Susceptible", "Exposed", "Infective", "Recovered"))))
  summ[1,1] = mean(out$HS1) + mean(out$HS2) + mean(out$HS3)
  summ[1,2] = mean(out$HE1) + mean(out$HE2) + mean(out$HE3)
  summ[1,3] = mean(out$HI1) + mean(out$HI2) + mean(out$HI3)
  summ[1,4] = mean(out$HR1) + mean(out$HR2) + mean(out$HR3)
  summ[2,3] = max(out$HI1 + out$HI2 + out$HI3)
  summ[2,4] = mean(out$HR1+out$HR2+out$HR3)/mean(out$HS1+out$HS2+out$HS3+out$HE1+out$HE2+out$HE3+out$HI1+out$HI2+out$HI3+out$HR1+out$HR2+out$HR3)
  summ[3,1] = mean(out$MS1) + mean(out$MS2) + mean(out$MS3)
  summ[3,2] = mean(out$ME1) + mean(out$ME2) + mean(out$ME3)
  summ[3,3] = mean(out$MI1) + mean(out$MI2) + mean(out$MI3)
  summ[3,4] = mean(out$MR1) + mean(out$MR2) + mean(out$MR3)
  summ[4,3] = max(out$MI1 + out$MI2 + out$MI3)
  summ[4,4] = mean(out$MR1+out$MR2+out$MR3)/mean(out$MS1+out$MS2+out$MS3+out$ME1+out$ME2+out$ME3+out$MI1+out$MI2+out$MI3+out$MR1+out$MR2+out$MR3)
  output$summary = renderTable(summ, digits = 2, rownames = T, na = "")
  ### Vectors
  summ2 = as.data.frame(array(NA, c(4,4), dimnames = list(c("Vector A", "Vector B", "Vector C", "Vector D"), c("Clean eggs", "Infected eggs", "Susc. adult", "Inf. adult"))))
  summ2[1,1] = mean(out$AP1) + mean(out$AP2) + mean(out$AP3)
  summ2[1,2] = mean(out$AQ1) + mean(out$AQ2) + mean(out$AQ3)
  summ2[1,3] = mean(out$AS1) + mean(out$AS2) + mean(out$AS3)
  summ2[1,4] = mean(out$AI1) + mean(out$AI2) + mean(out$AI3)
  summ2[2,1] = mean(out$BP1) + mean(out$BP2) + mean(out$BP3)
  summ2[2,2] = mean(out$BQ1) + mean(out$BQ2) + mean(out$BQ3)
  summ2[2,3] = mean(out$BS1) + mean(out$BS2) + mean(out$BS3)
  summ2[2,4] = mean(out$BI1) + mean(out$BI2) + mean(out$BI3)
  summ2[3,1] = mean(out$CP1) + mean(out$CP2) + mean(out$CP3)
  summ2[3,3] = mean(out$CS1) + mean(out$CS2) + mean(out$CS3)
  summ2[3,4] = mean(out$CI1) + mean(out$CI2) + mean(out$CI3)
  summ2[4,1] = mean(out$DP1) + mean(out$DP2) + mean(out$DP3)
  summ2[4,3] = mean(out$DS1) + mean(out$DS2) + mean(out$DS3)
  summ2[4,4] = mean(out$DI1) + mean(out$DI2) + mean(out$DI3)
  output$summary2 = renderTable(summ2, digits = 2, rownames = T, na = "")
  if(year>=7){
    summ3 = as.data.frame(array(NA, c(2,3), dimnames = list(c("Human", "Mammal"), c("Year +2", "Year +4", "Year +6"))))
    ss = ((year-7) %/% 10) * 36000; s3 = ss + 2*3600+1; e3 = ss + 3*3600; s5 = ss + 4*3600+1; e5 = ss + 5*3600; s7 = ss + 6*3600+1; e7 = ss + 7*3600
    summ3[1,1] = mean(out$HR1[s3:e3]+out$HR2[s3:e3]+out$HR3[s3:e3])/mean(out$HS1[s3:e3]+out$HS2[s3:e3]+out$HS3[s3:e3]+out$HE1[s3:e3]+out$HE2[s3:e3]+out$HE3[s3:e3]+out$HI1[s3:e3]+out$HI2[s3:e3]+out$HI3[s3:e3]+out$HR1[s3:e3]+out$HR2[s3:e3]+out$HR3[s3:e3])
    summ3[1,2] = mean(out$HR1[s5:e5]+out$HR2[s5:e5]+out$HR3[s5:e5])/mean(out$HS1[s5:e5]+out$HS2[s5:e5]+out$HS3[s5:e5]+out$HE1[s5:e5]+out$HE2[s5:e5]+out$HE3[s5:e5]+out$HI1[s5:e5]+out$HI2[s5:e5]+out$HI3[s5:e5]+out$HR1[s5:e5]+out$HR2[s5:e5]+out$HR3[s5:e5])
    summ3[1,3] = mean(out$HR1[s7:e7]+out$HR2[s7:e7]+out$HR3[s7:e7])/mean(out$HS1[s7:e7]+out$HS2[s7:e7]+out$HS3[s7:e7]+out$HE1[s7:e7]+out$HE2[s7:e7]+out$HE3[s7:e7]+out$HI1[s7:e7]+out$HI2[s7:e7]+out$HI3[s7:e7]+out$HR1[s7:e7]+out$HR2[s7:e7]+out$HR3[s7:e7])
    summ3[2,1] = mean(out$MR1[s3:e3]+out$MR2[s3:e3]+out$MR3[s3:e3])/mean(out$MS1[s3:e3]+out$MS2[s3:e3]+out$MS3[s3:e3]+out$ME1[s3:e3]+out$ME2[s3:e3]+out$ME3[s3:e3]+out$MI1[s3:e3]+out$MI2[s3:e3]+out$MI3[s3:e3]+out$MR1[s3:e3]+out$MR2[s3:e3]+out$MR3[s3:e3])
    summ3[2,2] = mean(out$MR1[s5:e5]+out$MR2[s5:e5]+out$MR3[s5:e5])/mean(out$MS1[s5:e5]+out$MS2[s5:e5]+out$MS3[s5:e5]+out$ME1[s5:e5]+out$ME2[s5:e5]+out$ME3[s5:e5]+out$MI1[s5:e5]+out$MI2[s5:e5]+out$MI3[s5:e5]+out$MR1[s5:e5]+out$MR2[s5:e5]+out$MR3[s5:e5])
    summ3[2,3] = mean(out$MR1[s7:e7]+out$MR2[s7:e7]+out$MR3[s7:e7])/mean(out$MS1[s7:e7]+out$MS2[s7:e7]+out$MS3[s7:e7]+out$ME1[s7:e7]+out$ME2[s7:e7]+out$ME3[s7:e7]+out$MI1[s7:e7]+out$MI2[s7:e7]+out$MI3[s7:e7]+out$MR1[s7:e7]+out$MR2[s7:e7]+out$MR3[s7:e7])
    output$summary3 = renderTable(summ3, digits = 3, rownames = T, na = "")
  }
} # if input$runIt
}) # observe
}
options(warn=0)

drawPlot4 = function(tt, g1, g2, g3, g4, y_ax1, y_ax2, legText, ps, pe, year)
{
  ps = (ps - 1)*0.1; pe = (pe - 1)*0.1
  if(pe%%30==0) pe = pe - 1
  ll = rep(c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"),
           (year+1))[(ps%/%30+1):(pe%/%30+1)]
  par(mar=c(5,5,1,5)); ymax1 = max(g1, g2); ymax2 = max(g3, g4)
  plot(tt, g1, type="l", col="blue", ylim=c(0,ymax1), ylab=y_ax1, xlab="Days", xaxt = "n")
  lines(tt, g2, col = "green")
  par(new=TRUE)
  plot(tt, g3, type="l", col="orange", ylim=c(0,ymax2), axes = F, ylab = NA, xlab=NA)
  axis(4); mtext(y_ax2,side=4,line=3)
  axis(1, at = c(-15+30*((ps%/%30+1):(pe%/%30+1))), labels = ll)
  lines(tt, g4, col="red")
  legend("topright", legText, lty=c(1,1,1,1), lwd=c(2,2,2,2), col=c("blue", "green", "orange", "red"))
}

drawPlot3 = function(tt, g1, g2, g3, y_ax1, y_ax2, legText, ps, pe, year)
{
  ps = (ps - 1)*0.1; pe = (pe - 1)*0.1
  if(pe%%30==0) pe = pe - 1
  ll = rep(c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"),
           (year+1))[(ps%/%30+1):(pe%/%30+1)]
  par(mar=c(5,5,1,5)); ymax1 = max(g1, g2); ymax2 = max(g3)
  plot(tt, g1, type="l", col="blue", ylim=c(0,ymax1), ylab=y_ax1, xlab="Days", xaxt = "n")
  lines(tt, g2, col = "green")
  par(new=TRUE)
  plot(tt, g3, type="l", col="red", ylim=c(0,ymax2), axes = F, ylab = NA, xlab=NA)
  axis(4); mtext(y_ax2,side=4,line=3)
  axis(1, at = c(-15+30*((ps%/%30+1):(pe%/%30+1))), labels = ll)
  legend("topright", legText, lty=c(1,1,1,1), lwd=c(2,2,2,2), col=c("blue", "green", "red"))
}

shinyApp(ui = ui, server = server)