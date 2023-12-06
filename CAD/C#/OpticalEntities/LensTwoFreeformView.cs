using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Controls;
using MathKernel;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    public class LensTwoFreeformView : SolidOpticalEntityView
    {
        // Field
        private LensTwoFreeform modelObject
        {
            get { return (LensTwoFreeform)base.ModelObject; }
        }
        private LensTwoFreeformZAView zaParams;
        private double skirtHeight;
        public static int numSurfaces = 5;

        // Properties
        public LensTwoFreeformZAView ZAParams
        {
            get
            {
                if (zaParams == null)
                    zaParams = new LensTwoFreeformZAView(modelObject.ZAParams, Material);
                return zaParams;
            }
        }
        public double SkirtHeight
        {
            get { return skirtHeight; }
            set { SetProperty<LensTwoFreeformView, double>(ref skirtHeight, value); }
        }
        public int NumSurfaces
        {
            get { return numSurfaces; }
        }

        public LensTwoFreeformView(LensTwoFreeform lens)
            :base(lens)
        {
            skirtHeight = lens.SkirtHeight;

            SubscribeToModelObject<LensTwoFreeform, LensTwoFreeformView>(lens);
        }

        static LensTwoFreeformView()
        {
            InitRules();
        }

        private static void InitRules()
        {
            RuleBook<LensTwoFreeformView>.ValidationRules.Add(new DelegateValidationRule<LensTwoFreeformView>(
                "SkirtHeight",
                "SkirtHeightOfRange",
                x => Epsilon.CompareNumerics(x.SkirtHeight, 0) > 0)
                );
            RuleBook<LensTwoFreeformView>.BindRules.Add(new PropertyBindRule<LensTwoFreeformView>
                (
                "ZAParams",
                "ZAParams",
                null,
                x => x.zaParams = null
            ));
            RuleBook<LensTwoFreeformView>.BindRules.Add(new PropertyBindRule<LensTwoFreeformView>
                (
                "SkirtHeight",
                "SkirtHeight",
                x => x.modelObject.SkirtHeight =x.SkirtHeight,
                x => x.skirtHeight = x.modelObject.SkirtHeight
            ));
        }
    }
}
