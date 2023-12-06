using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathKernel;
using ObjectModel.LightDistributions;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.Light_Distributions;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    using RuleBook = RuleBook<LensAxisymTIRFreeformView>;
    using ValidationRule = DelegateValidationRule<LensAxisymTIRFreeformView>;
    using BindRule = PropertyBindRule<LensAxisymTIRFreeformView>;

    public class LensAxisymTIRFreeformView : SolidOpticalEntityView
    {
        // Field
        private LensAxisymTIRFreeform modelObject;

        private LensAxisymTIRFreeformZAView zaParams;
        private static int numSurfaces = 6;
        private double neckHeight;
        // Properties
        public LensAxisymTIRFreeformZAView ZAParams
        {
            get
            {
                if (zaParams == null) zaParams = new LensAxisymTIRFreeformZAView(modelObject.ZAParams, Material);
                return zaParams;
            }
        }
        public double NeckHeight
        {
            get { return neckHeight; }
            set { SetProperty<LensAxisymTIRAsphericView, double>(ref neckHeight, value); }
        }

        public int NumSurfaces
        {
            get { return numSurfaces; }
        }


        public LensAxisymTIRFreeformView(LensAxisymTIRFreeform lens)
            : base(lens)
        {
            neckHeight = lens.NeckHeight;

            SubscribeToModelObject<LensAxisymTIRFreeform, LensAxisymTIRFreeformView>(lens);
        }

        static LensAxisymTIRFreeformView()
        {
            InitRules();
        }

        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
            "NeckHeight",
            "NeckHeightOfRange",
            x => Epsilon.CompareNumerics(x.NeckHeight, 0) > 0)
            );

            RuleBook.BindRules.Add(
                new BindRule(
                "ZAParams",
                "ZAParams",
                null,
                x => x.zaParams = null
                ));

            RuleBook.BindRules.Add(
                new BindRule(
                "NeckHeight",
                "NeckHeight",
                x => x.modelObject.NeckHeight = x.NeckHeight,
                x => x.neckHeight = x.modelObject.NeckHeight
                ));
        }
    }
}
