using System;
using MathKernel;
using ObjectModel.OpticalSceneEntities;
using ViewModelCore.Common;

namespace ViewModelCore.OpticalSceneEntities
{
    using RuleBook = RuleBook<EmittableSolidView>;
    using ValidationRule = DelegateValidationRule<EmittableSolidView>;
    using BindRule = PropertyBindRule<EmittableSolidView>;

    public enum EmittableSolidType
    {
        Cylinder,
        Parallelepiped
    };

    public abstract class EmittableSolidView : NotifyDataErrorInfo
    {
        internal abstract EmittableSolid ModelObject { get; }

        public abstract EmittableSolidType Type { get; }

        internal static EmittableSolidView GetSolidByType(EmittableSolidType type)
        {
            switch (type)
            {
                case EmittableSolidType.Cylinder:
                    return new CylinderView(new Cylinder());
                case EmittableSolidType.Parallelepiped:
                    return new ParallelepipedView(new Parallelepiped());
                default:
                    return null;
            }
        }

        internal static EmittableSolidView Wrap(EmittableSolid solid)
        {
            if (solid is Parallelepiped)
            {
                return new ParallelepipedView((Parallelepiped) solid);
            }
            if (solid is Cylinder)
            {
                return new CylinderView((Cylinder) solid);
            }

            throw new NotImplementedException();
        }
    }
}
