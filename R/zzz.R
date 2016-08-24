.onAttach <- function(...) {
  packageStartupMessage("\nUse 'statTarget.gui()' to restart the programe.\n",fill=TRUE)
  #statTarget::statTargetGUI()
  if (interactive()) statTarget.gui()
}
