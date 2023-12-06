namespace ViewModelCore.OpticalSceneEntities
{
    internal interface IEmittableSolidContainer
    {
        EmittableSolidType SolidType { get; set; }

        EmittableSolidView Solid { get; set; }
    }
}
