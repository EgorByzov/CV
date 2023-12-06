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
    using RuleBook = RuleBook<LensAxisymTIRAsphericView>;
    using ValidationRule = DelegateValidationRule<LensAxisymTIRAsphericView>;
    using BindRule = PropertyBindRule<LensAxisymTIRAsphericView>;
    public class LensAxisymTIRAsphericView : SolidOpticalEntityView
    {
        // Field
        private LensAxisymTIRAspheric modelObject
        {
            get { return (LensAxisymTIRAspheric)base.ModelObject; }
        }

        private double neckHeight;
        private LensAxisymTIRAsphericZAView zaParams;
        public static int numSurfaces = 6;

        // Properties
        public LensAxisymTIRAsphericZAView ZAParams
        {
            get
            {
                if (zaParams == null) zaParams = new LensAxisymTIRAsphericZAView(modelObject.ZAParams, Material);
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
        //constructors
        public LensAxisymTIRAsphericView(LensAxisymTIRAspheric lens)
            :base(lens)
        {
            neckHeight = lens.NeckHeight;

            SubscribeToModelObject<LensAxisymTIRAspheric, LensAxisymTIRAsphericView>(lens);
        }

        static LensAxisymTIRAsphericView()
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
