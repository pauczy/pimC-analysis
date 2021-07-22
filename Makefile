
APP_NAME      := analysis
SOURCE_FILES  := main.cc analysis.cc 
INSTALL_DIR   := .

USES_RFIO     := no
USES_ORACLE   := yes
USES_GFORTRAN := yes
HOME:=/lustre/hades/user/shower/Paulina_Czyz/analysis


include $(HADDIR)/hades.def.mk

#LIB_DIRS += $(PLUTODIR) ${HOME}/usr/lib64
#INC_DIRS += ${HOME}/usr/include
HYDRA_LIBS    += -lDst


.PHONY:  default
#default: clean build install
default: build install

include $(HADDIR)/hades.app.mk

