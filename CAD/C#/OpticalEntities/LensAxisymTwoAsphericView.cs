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
    using RuleBook = RuleBook<LensAxisymTwoAsphericView>;
    using ValidationRule = DelegateValidationRule<LensAxisymTwoAsphericView>;
    using BindRule = PropertyBindRule<LensAxisymTwoAsphericView>;

    public class LensAxisymTwoAsphericView : SolidOpticalEntityView
    {
        // Field
        private LensAxisymTwoAspheric modelObject
        {
            get { return (LensAxisymTwoAspheric) base.ModelObject; }
        }
        private LensAxisymTwoAsphericZAView zaParams;
        public static int numSurfaces = 3;
        private double skirtHeight;

        // Properties
        public LensAxisymTwoAsphericZAView ZAParams
        {
            get
            {
                if (zaParams == null) zaParams = new LensAxisymTwoAsphericZAView(modelObject.ZAParams, Material);
                return zaParams;
            }
        }

        public double SkirtHeight
        {
            get { return skirtHeight; }
            set { SetProperty<LensAxisymTwoAsphericView, double>(ref skirtHeight, value); }
        }
        public int NumSurfaces
        {
            get { return numSurfaces; }
        }

        public LensAxisymTwoAsphericView(LensAxisymTwoAspheric lens)
            :base(lens)
        {
            skirtHeight = lens.SkirtHeight;

            SubscribeToModelObject<LensAxisymTwoAspheric, LensAxisymTwoAsphericView>(lens);
        }

        static LensAxisymTwoAsphericView()
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
