if(PERIDYNO_LIBRARY_STRESSAWARESPH AND PERIDYNO_QT_GUI)
    set(LIB_DEPENDENCY     
        Core 
        Framework 
        IO          
        QtGUI
        StressAwareSPH
        SemiAnalyticalScheme
        ParticleSystem
        GLRenderEngine)
    add_example(Qt_StressAwareSPH_FishAndBall StressAwareSPH LIB_DEPENDENCY)
endif()