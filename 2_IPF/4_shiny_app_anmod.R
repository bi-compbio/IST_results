library(shiny)

# launch app animal models
dir_current <- getwd()
shinyOptions(ist.result.path = paste0(dir_current, "/2_IPF/3_IST_disease_output/ist_res_anmod.rds"))
shinyOptions(ist.browser.title = "IPF - animal models")
devtools::load_all("../istbrowser")
shiny::shinyAppDir(system.file("app", package = "ISTBrowser"))
