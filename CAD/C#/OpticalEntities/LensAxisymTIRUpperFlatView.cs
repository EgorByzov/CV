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
    using RuleBook = RuleBook<LensAxisymTIRUpperFlatView>;
    using ValidationRule = DelegateValidationRule<LensAxisymTIRUpperFlatView>;
    using BindRule = PropertyBindRule<LensAxisymTIRUpperFlatView>;

    public class LensAxisymTIRUpperFlatView : SolidOpticalEntityView
    {
        // Field
        private LensAxisymTIRUpperFlat modelObject
        {
            get { return (LensAxisymTIRUpperFlat)base.ModelObject; }
        }

        private double neckHeight;
        private LensAxisymTIRUpperFlatZAView zaParams;
        public static int numSurfaces = 5;

        // Properties
        public LensAxisymTIRUpperFlatZAView ZAParams
        {
            get
            {
                if (zaParams == null) zaParams = new LensAxisymTIRUpperFlatZAView(modelObject.ZAParams, Material);
                return zaParams;
            }
        }
        public double NeckHeight
        {
            get { return neckHeight; }
            set { SetProperty<LensAxisymTIRUpperFlatView, double>(ref neckHeight, value); }
        }
        public int NumSurfaces
        {
            get { return numSurfaces; }
        }

        public LensAxisymTIRUpperFlatView(LensAxisymTIRUpperFlat lens)
            :base(lens)
        {
            neckHeight = lens.NeckHeight;
            SubscribeToModelObject<LensAxisymTIRUpperFlat, LensAxisymTIRUpperFlatView>(lens);
        }

        static LensAxisymTIRUpperFlatView()
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
