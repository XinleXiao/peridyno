if(PERIDYNO_LIBRARY_STRESSAWARESPH)
    set(LIB_DEPENDENCY     
        Core 
        Framework 
        IO          
        StressAwareSPH
        GlfwGUI
        GLRenderEngine)
    add_example(StressAwareSPH_Fountain StressAwareSPH LIB_DEPENDENCY)
endif()