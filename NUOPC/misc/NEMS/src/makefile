SHELL = /bin/sh -x

# variables that need to move between makefiles on different levels of the build tree
# via shell environment variables
export TOP := $(shell pwd)
export PHYS_INC

# check if configure.nems file exists

exist := $(wildcard ./conf/configure.nems)
ifneq ($(strip $(exist)),)
  include ./conf/configure.nems
else
  $(error Error: !!!!!! Configure file is missing. Please create ./conf/configure.nems  by running ./configure !!!!!!)
endif

# conditionally read in the paths for external components

ifndef NEMSAppBuilder
  exist := $(wildcard ./conf/externals.nems)
  ifneq ($(strip $(exist)),)
    include ./conf/externals.nems
  else
    $(error Error: !!!!!! Externals file is missing. Please create ./conf/externals.nems  by running ./configure !!!!!!)
  endif
endif

# Need a variable that contains only a comma for some of the operations below
COMMA := ,

# external COMP dependencies ###################################################
# ATMs
ifneq (,$(findstring $(COMMA)satm$(COMMA),$(COMP)))
  include $(satm_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_SATM=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)xatm$(COMMA),$(COMP)))
  include $(xatm_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_XATM=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)datawam$(COMMA),$(COMP)))
  include $(datawam_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_DATAWAM=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
# OCNs
ifneq (,$(findstring $(COMMA)socn$(COMMA),$(COMP)))
  include $(socn_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_SOCN=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)xocn$(COMMA),$(COMP)))
  include $(xocn_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_XOCN=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)hycom$(COMMA),$(COMP)))
  include $(hycom_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_HYCOM=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)mom5$(COMMA),$(COMP)))
  include $(mom5_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_MOM5=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)pom$(COMMA),$(COMP)))
  include $(pom_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_POM=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
# SEAICEs
ifneq (,$(findstring $(COMMA)sice$(COMMA),$(COMP)))
  include $(sice_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_SICE=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)xice$(COMMA),$(COMP)))
  include $(xice_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_XICE=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)cice$(COMMA),$(COMP)))
  include $(cice_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_CICE=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)sis2$(COMMA),$(COMP)))
  include $(sis2_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_SIS2=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
# WAVs
ifneq (,$(findstring $(COMMA)swav$(COMMA),$(COMP)))
  include $(swav_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_SWAV=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)xwav$(COMMA),$(COMP)))
  include $(xwav_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_XWAV=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)ww3$(COMMA),$(COMP)))
  include $(ww3_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_WW3=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
# LNDs
ifneq (,$(findstring $(COMMA)slnd$(COMMA),$(COMP)))
  include $(slnd_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_SLND=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)xlnd$(COMMA),$(COMP)))
  include $(xlnd_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_XLND=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)noah$(COMMA),$(COMP)))
  include $(noah_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_NOAH=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)noahmp$(COMMA),$(COMP)))
  include $(noahmp_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_NOAHMP=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)lis$(COMMA),$(COMP)))
  include $(lis_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_LIS=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
# IPMs
ifneq (,$(findstring $(COMMA)sipm$(COMMA),$(COMP)))
  include $(sipm_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_SIPM=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)xipm$(COMMA),$(COMP)))
  include $(xipm_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_XIPM=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)ipe$(COMMA),$(COMP)))
  include $(ipe_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_IPE=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)dataipe$(COMMA),$(COMP)))
  include $(dataipe_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_DATAIPE=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
# HYDs
ifneq (,$(findstring $(COMMA)shyd$(COMMA),$(COMP)))
  include $(shyd_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_SHYD=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)xhyd$(COMMA),$(COMP)))
  include $(xhyd_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_XHYD=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifneq (,$(findstring $(COMMA)wrfhydro$(COMMA),$(COMP)))
  include $(wrfhydro_mk)
  DEP_FRONTS    := $(DEP_FRONTS) -DFRONT_WRFHYDRO=$(ESMF_DEP_FRONT)
  DEP_INCS      := $(DEP_INCS) $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS := $(DEP_CMPL_OBJS) $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS := $(DEP_LINK_OBJS) $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH := $(DEP_SHRD_PATH) $(addprefix -L, $(ESMF_DEP_SHRD_PATH)) $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS := $(DEP_SHRD_LIBS) $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
################################################################################

# development with NEMS phys to be compiled every time. Assumes std phys location
# no pre-build libphys.a will be used 
#  PHYS_MODE must be set as "compile"

ifeq ($(PHYS_MODE),compile)
	PHYS_LIB = $(TOP)/atmos/gsm/gsmphys
	PHYS_INC = $(TOP)/atmos/gsm/gsmphys
	PHYS_DIR = $(TOP)/atmos/gsm/gsmphys
endif
ifeq ($(CHEM_MODE),compile)
        CHEM_LIB = $(TOP)/chem
        CHEM_INC = $(TOP)/chem/gocart/src/Config/
        CHEM_DIR = $(TOP)/chem
	CHEM_MOD = $(TOP)/chem/gocart/${ARCH}/include
ESMADIR=chem/gocart
include $(ESMADIR)/src/Config/NCEP_base.mk
endif


CPPFLAGS += $(DEP_FRONTS)

TARGET = ../exe/NEMS.x

INCS = -Iatmos \
	-I$(CHEM_MOD) \
       -Iatmos/share \
       -IENS_Cpl \
       -Iatmos/nmm \
       -Iatmos/fim \
       -I$(gfsdir) \
       -I$(PHYS_INC) \
	-I. \
       -Iatmos/gen \
	-Iatmos/share \
       $(DEP_INCS)

MAIN = MAIN_NEMS.o

OBJS = module_MEDIATOR_methods.o \
       module_MEDIATOR.o \
       module_MEDIATOR_SpaceWeather.o \
       module_EARTH_INTERNAL_STATE.o \
       module_EARTH_GRID_COMP.o \
       module_NEMS_INTERNAL_STATE.o \
       module_NEMS_GRID_COMP.o

COBJS = # c_code.o

LIBS = atmos/libatmos.a                 \
       atmos/nmm/libnmm.a               \
       $(gfsdir)/libgfs.a               \
       atmos/io/libio.a                 \
       atmos/gen/libgen.a               \
       atmos/fim/libfimtop.a            \
       $(gfsdir)/libutil/gfs_libutil.a  \
       $(gfsdir)/phys/gfs_physics.a     \
       $(PHYS_LIB)/libphys.a             \
       $(gfsdir)/dyn/gfs_dynamics.a     \
       atmos/share/libshare.a           \
       atmos/post/libpost.a             \
       $(DEP_LINK_OBJS)                 \
       ENS_Cpl/ENS_Cpl.a		\
       $(DEP_SHRD_PATH)	                \
       $(DEP_SHRD_LIBS)

LIBS_NMM = atmos/libatmos.a             \
           atmos/nmm/libnmm.a           \
           $(gfsdir)/libgfs.a           \
           atmos/fim/libfimtop.a        \
           atmos/gen/libgen.a           \
           $(PHYS_LIB)/libphys.a         \
           atmos/share/libshare.a       \
           atmos/post/libpost.a         \
           $(DEP_LINK_OBJS)             \
           ENS_Cpl/ENS_Cpl.a

LIBS_GEN = atmos/libatmos.a             \
           atmos/nmm/libnmm.a           \
           atmos/io/libio.a             \
           $(gfsdir)/libgfs.a           \
           atmos/fim/libfimtop.a        \
           atmos/gen/libgen.a           \
           $(PHYS_LIB)/libphys.a         \
           atmos/share/libshare.a       \
           atmos/post/libpost.a         \
           $(DEP_LINK_OBJS)             \
           ENS_Cpl/ENS_Cpl.a

#TODO:  DRY these settings, they are redundant with atmos/fim/makefile and 
#TODO:  atmos/fim/FIMsrc ...  
SMS      = /gpfs/t3/global/save/wx20tbh/SMS/SMS_r237
SMSLIB   = $(SMS)/lib/libsms.a
# libraries created during build inside atmos/fim/FIM_src_devccs
FIMSRCLIBS = atmos/fim/FIMsrc_devccs/lib/libfim.a      \
             atmos/fim/FIMsrc_devccs/lib/libwrfio.a    \
             atmos/fim/FIMsrc_devccs/lib/libgribio.a   \
             atmos/fim/FIMsrc_devccs/lib/libvlint.a    \
             atmos/fim/FIMsrc_devccs/lib/libsp_4.a     \
             atmos/fim/FIMsrc_devccs/lib/libslint.a    \
             atmos/fim/FIMsrc_devccs/lib/libsigio_4.a  \
             atmos/fim/FIMsrc_devccs/lib/libsfcio_4.a  \
             atmos/fim/FIMsrc_devccs/lib/libw3_4.a     \
             atmos/fim/FIMsrc_devccs/lib/libbacio_4.a  \
             atmos/fim/FIMsrc_devccs/lib/libcntl.a     \
             atmos/fim/FIMsrc_devccs/lib/libsysshare.a

LIBS_FIM = atmos/libatmos.a                     \
           atmos/nmm/libnmm.a                   \
           $(gfsdir)/libgfs.a                   \
           atmos/fim/libfimtop.a                \
           $(SMSLIB) $(FIMSRCLIBS)              \
           atmos/gen/libgen.a                   \
           atmos/share/libshare.a               \
           $(DEP_LINK_OBJS)                     \
           ENS_Cpl/ENS_Cpl.a                    \
           -L/nwprod/lib -lnetcdf_64

MAKEFILE = makefile

#
# GOCART specific compilation variables/flags
#
GOCART_MODE=stub
ESMADIR=$(CHEM_DIR)
include $(CHEM_INC)/NCEP_base.mk

ifeq ($(GOCART_MODE),stub)
	LIBGOCART=
else
	LIBGOCART=$(LIB_GOCART)
		CHEM_LIB=$(CHEM_DIR)_full/$(CHEM_VER)
		CHEM_MOD=$(CHEM_LIB)/mod
		CHEM_INC=$(CHEM_LIB)/include
		PHYS_LIB=$(PHYS_DIR)_full/$(CHEM_VER)
		PHYS_MOD=$(PHYS_LIB)/mod
		PHYS_INC=$(PHYS_LIB)/include
endif

ifeq ($(PHYS_MODE),compile)
        PHYS_LIB = $(TOP)/atmos/gsm/gsmphys
        PHYS_INC = $(TOP)/atmos/gsm/gsmphys
        PHYS_DIR = $(TOP)/atmos/gsm/gsmphys
endif
ifeq ($(CHEM_MODE),compile)
        CHEM_LIB = $(TOP)/chem/gocart
        CHEM_INC = $(TOP)/chem/gocart/src/Config
        CHEM_DIR = $(TOP)/chem
	CHEM_MOD = $(TOP)/chem/gocart/${ARCH}/include
ESMADIR=chem/gocart
include $(ESMADIR)/src/Config/NCEP_base.mk
endif



export GOCART_MODE
ESMADIR=$(CHEM_LIB)
include $(CHEM_INC)/NCEP_base.mk
#
# log info
#
WHOAMI  := $(shell whoami)
CMP_DATE := $(shell date)
SVN_REV := $(shell svnversion -n .)
SVN_URL := $(shell svn info | grep URL | cut -f 2 -d ' ')
CMP_YEAR := $(shell date +"%Y" )
CMP_JD := $(shell date +"%j" )

CPPFLAGS += -D'SVN_INFO="($(WHOAMI)) $(CMP_DATE) r$(SVN_REV) $(SVN_URL)"'

ifdef CMP_YEAR
  CPPFLAGS += -D'CMP_YEAR=$(CMP_YEAR)'
endif
ifdef CMP_JD
  CPPFLAGS += -D'CMP_JD=$(CMP_JD)'
endif

################################################################################
#
# Ask user which target he/she wants to build
#
################################################################################
all:
	@echo "       CHEM_MODi----------->" $(CHEM_MOD)
	@echo " gmake nmm_gsm ----------------------------- NMM and GSM - without the GOCART"
	@echo " gmake gen --------------------------------- GEN "
	@echo " gmake gsm_gen_post GOCART_MODE=full ------- GEN and GSM - with the GOCART and post"
	@echo " gmake nmm_gsm GOCART_MODE=full ------------ NMM and GSM - with the GOCART"
	@echo " gmake nmm --------------------------------- NMM "
	@echo " gmake nmm_post ---------------------------- NMM - with post"
	@echo " gmake gsm --------------------------------- GSM - without the GOCART"
	@echo " gmake gsm_post ---------------------------- GSM - with post without the GOCART"
	@echo " gmake gsm GOCART_MODE=full ---------------- GSM - with the GOCART"
	@echo " gmake fim --------------------------------- FIM "
#	@echo " gmake fim_gsm ----------------------------- FIM and GSM - without the GOCART"
#	@echo " gmake fim_gsm GOCART_MODE=full ------------ FIM and GSM - with the GOCART"
	@echo
	@echo "In order to build NEMS with external components, all of the above targets support"
	@echo "the optional COMP variable. Currently supported values are: 
	@echo "\"satm\", \"xatm\";"
	@echo "\"socn\", \"xocn\", \"mom5\", \"hycom\", \"pom\";"
	@echo "\"sice\", \"xice\", \"cice\";"
	@echo "\"swav\", \"xwav\", \"ww3\";"
	@echo "\"slnd\", \"xlnd\", \"noah\", \"noahmp\";"
	@echo "\"sipm\", \"xipm\", \"ipe\";"
	@echo "\"shyd\", \"xhyd\", \"wrfhydro\"."
	@echo
	@echo "E.g. \"gmake nmm_gsm COMP=satm,xatm,socn,xocn,mom5,hycom,sice,xice,cice,swav,xwav,slnd,xlnd,sipm,xipm,shyd,xhyd\"."
	@echo

################################################################################
#
# NMM
#
################################################################################

nmm: nmm_libs $(OBJS) $(COBJS) $(MAIN)
	$(FC) $(FFLAGS_NMM) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(COBJS) $(LIBS_NMM) $(EXTLIBS) 
	@echo "$(TARGET) is created for NMM core."

nmm_libs:
	@echo "       CHEM_MODi----------->" $(CHEM_MOD)
	cd atmos/share         && gmake $(J)
ifeq ($(PHYS_MODE),compile)
	cd atmos/gsm/gsmphys          && gmake $(J) gocart_stub
endif
	cd atmos/nmm           && gmake $(J)
	cd $(gfsdir)           && gmake stub
	cd atmos/gen           && gmake stub
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake $(J)
	cd ENS_Cpl             && gmake stub

################################################################################
#
# NMM with post
#
################################################################################

nmm_post: nmm_post_libs $(OBJS) $(COBJS) $(MAIN)
	$(FC) $(FFLAGS_NMM) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(COBJS) $(LIBS_NMM) $(EXTLIBS_POST) 
	@echo "$(TARGET) is created for NMM core."

nmm_post_libs:
	cd atmos/share         && gmake $(J)
ifeq ($(PHYS_MODE),compile)
	cd atmos/gsm/gsmphys          && gmake $(J) gocart_stub
endif
	cd atmos/nmm           && gmake $(J)
	cd $(gfsdir)           && gmake stub
	cd atmos/gen           && gmake stub
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake $(J) gfs_stub
	cd atmos               && gmake $(J)
	cd ENS_Cpl             && gmake stub

################################################################################
#
# GSM
#
################################################################################

gsm: gsm_libs $(OBJS) $(COBJS) $(MAIN)
	@echo "       CHEM_MODE ----------->"$(CHEM_MODE)
	@echo "       PHYS_MODE ----------->"$(PHYS_MODE)
	@echo "       PHYS_INC----------->"$(PHYS_INC)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(COBJS) $(LIBS) $(LIBGOCART) $(gfsdir)/phys/gfs_physics.a  $(gfsdir)/dyn/gfs_dynamics.a $(EXTLIBS)
	@echo "NEMS.x successfully created for GSM core."

gsm_libs:
ifeq ($(GOCART_MODE),full)
ifeq ($(CHEM_MODE),compile)
	cd chem/gocart/src     && gmake install GOCART_MODE=$(GOCART_MODE)
endif
endif
	@echo "       CHEM_MODE ----------->" $(CHEM_MODE)
	@echo "       PHYS_MODE ----------->" $(PHYS_MODE)
	 @echo "       PHYS_INC----------->"$(PHYS_INC)
	cd atmos/share         && gmake
ifeq ($(PHYS_MODE),compile)
	@echo "       CHEM_MODi----before phys------->" $(CHEM_MOD)
	cd atmos/gsm/gsmphys          && gmake gocart_$(GOCART_MODE)
endif
	cd atmos/nmm           && gmake stub
	cd atmos/gen           && gmake stub
	cd atmos/io            && gmake
	echo "       CHEM_MODi----------->ediri" $(CHEM_MOD)
	cd $(gfsdir)           && gmake gocart_$(GOCART_MODE)
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# GSM with post
#
################################################################################

gsm_post: gsm_post_libs $(OBJS) $(COBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(COBJS) $(LIBS) $(LIBGOCART) $(gfsdir)/phys/gfs_physics.a $(gfsdir)/dyn/gfs_dynamics.a $(EXTLIBS_POST)
	@echo "NEMS.x successfully created for GSM core."

gsm_post_libs:
ifeq ($(GOCART_MODE),full)
	ifeq ($(CHEM_MODE),compile)
	cd chem/gocart/src     && gmake install GOCART_MODE=$(GOCART_MODE)
	endif
endif
	cd atmos/share         && gmake
ifeq ($(PHYS_MODE),compile)
	cd atmos/gsm/gsmphys          && gmake gocart_$(GOCART_MODE)
endif
	cd atmos/nmm           && gmake stub
	cd atmos/gen           && gmake stub
	cd atmos/io            && gmake
	cd $(gfsdir)           && gmake gocart_$(GOCART_MODE)
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake nmm_stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# GEN
#
################################################################################

gen: gen_libs $(OBJS) $(COBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(COBJS) $(OBJS) $(LIBS_GEN) $(EXTLIBS)
	@echo "NEMS.x successfully created for GEN core."

gen_libs:
	cd atmos/share         && gmake
ifeq ($(PHYS_MODE),compile)
	cd atmos/gsm/gsmphys          && gmake
endif
	cd atmos/nmm           && gmake stub
	cd $(gfsdir)     && gmake stub
	cd atmos/fim           && gmake stub
	cd atmos/gen           && gmake
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake


################################################################################
#
# GEN, NMM and GSM
#
################################################################################

nmm_gsm_gen: nmm_gsm_gen_libs $(OBJS) $(COBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(COBJS) $(LIBS) $(LIBGOCART) $(gfsdir)/phys/gfs_physics.a $(gfsdir)/dyn/gfs_dynamics.a $(EXTLIBS)
	@echo "NEMS.x successfully created for GEN, NMM and GSM core."

nmm_gsm_gen_libs:
ifeq ($(GOCART_MODE),full)
	ifeq ($(CHEM_MODE),compile)
	cd chem/gocart/src     && gmake install GOCART_MODE=$(GOCART_MODE)
	endif
endif
	cd atmos/share         && gmake
ifeq ($(PHYS_MODE),compile)
	cd atmos/gsm/gsmphys          && gmake gocart_$(GOCART_MODE)
endif
	cd atmos/nmm           && gmake
	cd atmos/io            && gmake
	cd $(gfsdir)           && gmake gocart_$(GOCART_MODE)
	cd atmos/gen           && gmake 
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# GEN and GSM with post
#
################################################################################

gsm_gen_post: gsm_gen_post_libs $(OBJS) $(COBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(COBJS) $(LIBS) $(LIBGOCART) $(gfsdir)/phys/gfs_physics.a $(gfsdir)/dyn/gfs_dynamics.a $(EXTLIBS_POST)
	@echo "NEMS.x successfully created for GEN and GSM core."

gsm_gen_post_libs:
ifeq ($(GOCART_MODE),full)
	ifeq ($(CHEM_MODE),compile)
	cd chem/gocart/src     && gmake install GOCART_MODE=$(GOCART_MODE)
	endif
endif
	cd atmos/share         && gmake
ifeq ($(PHYS_MODE),compile)
	cd atmos/gsm/gsmphys          && gmake gocart_$(GOCART_MODE)
endif
	cd atmos/nmm           && gmake stub

	cd atmos/io            && gmake
	cd $(gfsdir)           && gmake gocart_$(GOCART_MODE)
	cd atmos/gen           && gmake
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake nmm_stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# NMM and GSM
#
################################################################################

nmm_gsm: nmm_gsm_libs $(OBJS) $(COBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(COBJS) $(OBJS) $(LIBS) $(LIBGOCART) $(gfsdir)/phys/gfs_physics.a $(gfsdir)/dyn/gfs_dynamics.a $(EXTLIBS)
	@echo "NEMS.x successfully created for NMM and GSM core."

nmm_gsm_libs:
ifeq ($(GOCART_MODE),full)
	ifeq ($(CHEM_MODE),compile)
	cd chem/gocart/src     && gmake install GOCART_MODE=$(GOCART_MODE)
	endif
endif
	cd atmos/share         && gmake
ifeq ($(PHYS_MODE),compile)
	cd atmos/gsm/gsmphys          && gmake gocart_$(GOCART_MODE)
endif
	cd atmos/nmm           && gmake
	cd atmos/io            && gmake
	cd $(gfsdir)           && gmake gocart_$(GOCART_MODE)
	cd atmos/fim           && gmake stub
	cd atmos/gen           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# FIM
#
################################################################################

# TODO:  Change to re-use NEMS standard $(EXTLIBS) instead of FIM versions 
# TODO:  currently in $(LIBS_FIM) .  
fim: fim_libs $(OBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS_FIM) $(ESMF_LIB)
	@echo "NEMS.x successfully created for FIM core."

fim_libs:
	cd atmos/share         && gmake
	cd atmos/nmm           && gmake stub
	cd $(gfsdir)     && gmake stub 
	cd atmos/fim           && gmake
	cd atmos/gen           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# FIM and GSM
#TBH:  not yet integrated with FIM build
#
################################################################################

fim_gsm: fim_gsm_libs $(OBJS) $(MAIN)
#	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS) $(EXTLIBS)
#	@echo "NEMS.x successfully created for FIM and GSM core."
	@echo "ERROR:  Cannot build FIM with other cores yet.  nems.x is NOT created for FIM and GSM core."
	@exit -1

fim_gsm_libs:
	cd atmos/share         && gmake
	cd atmos/nmm           && gmake stub
	cd atmos/io            && gmake
	cd $(gfsdir)           && gmake
	cd atmos/fim           && gmake
	cd atmos/gen           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# Clean everything
#
################################################################################

clean:
	cd atmos               && gmake clean
	cd ENS_Cpl             && gmake clean
ifeq ($(CHEM_MODE),compile)
		cd chem                && gmake clean GOCART_MODE=full
endif
	$(RM) -f *.f90 *.lst *.o *.mod lm map $(TARGET)

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

$(COBJS) : %.o: %.c 
	$(CC) -c $*.c

$(OBJS): %.o: %.f90
	$(FC) $(FFLAGS) $(INCS) -c $*.f90
	$(RM) $*.f90

MAIN_NEMS.o: %.o: %.f90
	$(FC) $(FFLAGS) $(INCS) -c $*.f90
	$(RM) $*.f90

