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
    using RuleBook = RuleBook<LensFreeformSphDomeView>;
    using ValidationRule = DelegateValidationRule<LensFreeformSphDomeView>;
    using BindRule = PropertyBindRule<LensFreeformSphDomeView>;

    public class LensFreeformSphDomeView : SolidOpticalEntityView
    {
        // Field
        private LensFreeformSphDome modelObject
        {
            get { return (LensFreeformSphDome)base.ModelObject; }
        }
        private LensFreeformSphDomeZAView zaParams;
        public static int numSurfaces = 5;
        private double skirtHeight;

        // Properties
        public LensFreeformSphDomeZAView ZAParams
        {
            get
            {
                if (zaParams == null) zaParams = new LensFreeformSphDomeZAView(modelObject.ZAParams, Material);
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

        public LensFreeformSphDomeView(LensFreeformSphDome lens)
            : base(lens)
        {
            skirtHeight = lens.SkirtHeight;

            SubscribeToModelObject<LensFreeformSphDome, LensFreeformSphDomeView>(lens);
        }

        static LensFreeformSphDomeView()
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
