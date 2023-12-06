using MathKernel;
using Newtonsoft.Json;
using ObjectModel.OpticalSceneEntities.OpticalEntities;
using ViewModelCore.Common;
using ViewModelCore.OpticalSceneEntities.OpticalEntities;

namespace ViewModelCore.OpticalSceneEntities
{

    using RuleBook = RuleBook<RegistratorView>;
    using ValidationRule = DelegateValidationRule<RegistratorView>;
    using BindRule = PropertyBindRule<RegistratorView>;

    [JsonObject(MemberSerialization.Fields)]
    public class RegistratorView : OpticalEntityView
    {
        //field
        public Registrator modelObject
        {
            get { return (Registrator) base.ModelObject; }
        }
        private double xSize;
        private double ySize;


        // properties
        public double XSize
        {
            get { return xSize; }
            set
            {
                SetProperty<RegistratorView, double>(ref xSize, value);
            }
        }
        public double YSize
        {
            get { return ySize; }
            set
            {
                SetProperty<RegistratorView, double>(ref ySize, value);
            }
        }

        //constructrs
        static RegistratorView()
        {
            InitRules();
        }
        
        [JsonConstructor]
        internal RegistratorView(Registrator modelObject)
            : base(modelObject)
        {
            xSize = modelObject.XSize;
            ySize = modelObject.YSize;
            SubscribeToModelObject<Registrator, RegistratorView>(modelObject);
        }

        private static void InitRules()
        {
            RuleBook.ValidationRules.Add(new ValidationRule(
                "XSize",
                "XSizeOutOfRange",
                x => Epsilon.CompareNumerics(x.XSize, 0) > 0)
                );
            RuleBook.ValidationRules.Add(new ValidationRule(
                "YSize",
                "YSizeOutOfRange",
                x => Epsilon.CompareNumerics(x.YSize, 0) > 0)
                );
            RuleBook.BindRules.Add(new BindRule(
                "XSize",
                "XSize",
                x => x.modelObject.XSize = x.XSize,
                x => x.xSize = x.modelObject.XSize
                ));
            RuleBook.BindRules.Add(new BindRule(
              "YSize",
              "YSize",
              x => x.modelObject.YSize = x.YSize,
              x => x.ySize = x.modelObject.YSize
              ));

        }
    }
}
