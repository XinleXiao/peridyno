if(PERIDYNO_LIBRARY_STRESSAWARESPH AND PERIDYNO_QT_GUI)
    set(LIB_DEPENDENCY     
        Core 
        Framework 
        IO          
        GlfwGUI
        StressAwareSPH
        GLRenderEngine)
    add_example(StressAwareSPH_CollidingFish StressAwareSPH LIB_DEPENDENCY)
endif()