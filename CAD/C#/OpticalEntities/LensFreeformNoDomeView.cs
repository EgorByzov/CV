using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathKernel;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    using RuleBook = RuleBook<LensFreeformNoDomeView>;
    using ValidationRule = DelegateValidationRule<LensFreeformNoDomeView>;
    using BindRule = PropertyBindRule<LensFreeformNoDomeView>;

    public class LensFreeformNoDomeView : SolidOpticalEntityView
    {
        // Field
        private LensFreeformNoDome modelObject
        {
            get { return (LensFreeformNoDome) base.ModelObject; }
        }
        private LensFreeformNoDomeZAView zaParams;
        public static int numSurfaces = 3;
        private double skirtHeight;
        // Properties
        public LensFreeformNoDomeZAView ZAParams
        {
            get
            {
                if (zaParams == null) zaParams = new LensFreeformNoDomeZAView(modelObject.ZAParams, Material);
                return zaParams;
            }
        }
        public double SkirtHeight
        {
            get { return skirtHeight; }
            set { SetProperty<LensFreeformSphDomeView, double>(ref skirtHeight, value); }
        }
        public int NumSurfaces
        {
            get { return numSurfaces; }
        }

        public LensFreeformNoDomeView(LensFreeformNoDome lens)
            :base(lens)
        {
            skirtHeight = lens.SkirtHeight;

            SubscribeToModelObject<LensFreeformNoDome, LensFreeformNoDomeView>(lens);
        }

        static LensFreeformNoDomeView()
        {
            InitRules();
        }

        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
            "SkirtHeight",
            "SkirtHeightOfRange",
            x => Epsilon.CompareNumerics(x.SkirtHeight, 0) > 0)
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
                "SkirtHeight",
                "SkirtHeight",
                x => x.modelObject.SkirtHeight = x.SkirtHeight,
                x => x.skirtHeight = x.modelObject.SkirtHeight
                ));
        }
    }
}
