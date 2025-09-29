source(here::here("deploy_mestools.R"))

check_mestools_status()


setup_upstream()


deploy_mestools_safe("Your commit message")
