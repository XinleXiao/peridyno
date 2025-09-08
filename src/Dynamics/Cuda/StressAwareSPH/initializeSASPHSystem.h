#pragma once
#include<Plugin/PluginEntry.h>

namespace dyno
{
	class SASPHInitializer :public PluginEntry
	{
	public:
		static PluginEntry* instance();

	protected:
		void initializeActions() override;

	private:
		SASPHInitializer();

		static std::atomic<SASPHInitializer*> gInstance;
		static std::mutex gMutex;
	};
}

namespace SASPH {
	dyno::PluginEntry* initStaticPlugin();

	PERIDYNO_API dyno::PluginEntry* initDynoPlugin();
}