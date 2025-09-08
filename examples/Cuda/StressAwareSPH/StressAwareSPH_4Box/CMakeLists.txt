if(PERIDYNO_LIBRARY_STRESSAWARESPH)
    set(LIB_DEPENDENCY     
        StressAwareSPH 
        Multiphysics
        GlfwGUI
        ABCExporter
    )
    add_example(StressAwareSPH_4Box StressAwareSPH LIB_DEPENDENCY)
endif()