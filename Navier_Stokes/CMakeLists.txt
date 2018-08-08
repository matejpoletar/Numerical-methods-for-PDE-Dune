add_executable("NS2" driver.hh bc_extension.hh exact.hh nrpdj-pdelab-evolNS.cc)
target_link_dune_default_libraries("NS2")
