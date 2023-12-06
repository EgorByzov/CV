using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.MaterialProperties;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    using RuleBook = RuleBook<SolidLensZAParamsView>;
    using ValidationRule = DelegateValidationRule<SolidLensZAParamsView>;
    using BindRule = PropertyBindRule<SolidLensZAParamsView>;

    public abstract class SolidLensZAParamsView : ZAParamsView, MaterialKeeper
    {
        //field
        private MaterialView innMaterial;
        private double wavelength;

        // Properties
        public abstract SolidLensZAParams ModelObject { get; }

        public MaterialView Material
        {
            get
            {
                if(innMaterial == null)
                    innMaterial = MaterialView.CreateByMaterial(ModelObject.InnMaterial);
                return innMaterial;
            }
            set
            {
                SetProperty(ref innMaterial, value);
            }
        }
        public double Wavelength
        {
            get { return wavelength; }
            set { SetProperty<SolidLensZAParamsView, double>(ref wavelength, value); }
        }

        //constructrs
        static SolidLensZAParamsView()
        {
            InitRules();
        }

        private static void InitRules()
        {   
            RuleBook.BindRules.Add(new BindRule(
                "Material",
                "InnMaterial",
                x =>
                {
                    x.ModelObject.InnMaterial = x.innMaterial.modelObject;
                },
                x => x.innMaterial = null
                ));

            RuleBook.BindRules.Add(new BindRule(
                "Wavelength",
                "Wavelength",
                x => x.ModelObject.Wavelength = x.Wavelength,
                x => x.Wavelength = x.ModelObject.Wavelength
                ));
        }

        public override abstract OpticalSceneObject ComputeOpticalElement(OpticalEntityView source,
            Light_Distributions.RequiredSimulationResultView reqSim);

    }
}
