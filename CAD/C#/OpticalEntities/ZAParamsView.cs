using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ObjectModel.LightDistributions;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.Light_Distributions;
using ViewModelCore.OpticalSceneEntities.OriginEntities;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    public abstract class ZAParamsView : NotifyDataErrorInfo
    {
        public abstract OpticalSceneObject ComputeOpticalElement(OpticalEntityView source,
            RequiredSimulationResultView reqSim);


    }
}
