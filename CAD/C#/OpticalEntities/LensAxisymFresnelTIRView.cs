using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;

namespace ViewModelCore.OpticalSceneEntities.OpticalEntities
{
    using RuleBook = RuleBook<LensAxisymFresnelTIRView>;
    using ValidationRule = DelegateValidationRule<LensAxisymFresnelTIRView>;
    using BindRule = PropertyBindRule<LensAxisymFresnelTIRView>;

    public class LensAxisymFresnelTIRView : SolidOpticalEntityView
    {
        // Field
        private LensAxisymFresnelTIR modelObject
        {
            get { return (LensAxisymFresnelTIR) base.ModelObject; }
        }
        private LensAxisymFresnelTIRZAView zaParams;
        private static int numSurfaces = 1;

        // Properties
        public LensAxisymFresnelTIRZAView ZAParams
        {
            get
            {
                if (zaParams == null) zaParams = new LensAxisymFresnelTIRZAView(modelObject.ZAParams, Material);
                return zaParams;
            }
        }
        public int NumSurfaces
        {
            get { return numSurfaces; }
        }
        public LensAxisymFresnelTIRView(LensAxisymFresnelTIR lens)
            :base(lens)
        {

            SubscribeToModelObject<LensAxisymFresnelTIR, LensAxisymFresnelTIRView>(lens);
        }

        static LensAxisymFresnelTIRView()
        {
            InitRules();
        }

        private static void InitRules()
        {
            RuleBook.BindRules.Add(
                new BindRule(
                "ZAParams",
                "ZAParams",
                null,
                x => x.zaParams = null
                ));
        }
    }
}
