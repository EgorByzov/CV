using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Collections.Specialized;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Media;
using System.Xml.Linq;
using Common;
using Common.TypeExtensions;
using ViewModelCore.Controllers;

namespace ViewModelCore.OpticalSceneEntities.SceneView
{
    internal enum ElementType { Registrator, Source }

    internal class OSOneLevelTreeView : TreeView, INotifyPropertyChanged
    {
        public event PropertyChangedEventHandler PropertyChanged;

        private Project _curProject;
        
        public ElementType ElementType { get; private set; }
        public OpenGLOpticalScene OSPresenter
        {
            get { return (_curProject == null) ? null : _curProject.OS.ScenePresenter; }
        }

        // ReSharper disable once MemberCanBeInternal
        public void InitCurrentProject(Project project, ElementType type)
        {
            Items.Clear();
            _curProject = project;

            if (_curProject == null) return;
        
            ObservableCollection<OpticalSceneObject> _elements =
                type == ElementType.Registrator ? _curProject.IlluminanceRegistrators : _curProject.Sources;
            
            foreach (var elem in _elements)
                Items.Add(new OSOneLevelTreeViewItem(elem, this));

            OnPropertyChanged("NumElements");

            _elements.CollectionChanged += (sender, args) =>
            {
                switch (args.Action)
                {
                    case NotifyCollectionChangedAction.Add:
                        Items.Add(new OSOneLevelTreeViewItem(args.NewItems[0] as OpticalSceneObject, this));
                        break;
                    case NotifyCollectionChangedAction.Remove:
                        foreach (var element in Items)
                            if (((OSOneLevelTreeViewItem) element).OSObject == args.OldItems[0])
                            {
                                Items.Remove(element);
                                break;
                            }
                        break;
                }
                OnPropertyChanged("NumElements");
            };
        }

        protected virtual void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            var handler = PropertyChanged;
            if (handler != null) handler(this, new PropertyChangedEventArgs(propertyName));
        }

        internal void UpdateColors()
        {
            foreach (OSOneLevelTreeViewItem item in Items)
                item.UpdateForeground();
        }
    }

    internal class OSOneLevelTreeViewItem : TreeViewItem
    {
        private OpticalSceneObject _osObject;
        private OSOneLevelTreeView _tree;

        public OpticalSceneObject OSObject {get { return _osObject; }}

        public OSOneLevelTreeViewItem(){}

        internal OSOneLevelTreeViewItem(OpticalSceneObject oso, OSOneLevelTreeView tree)
        {
            _osObject = oso;
            _tree = tree;

            SetBinding(HeaderProperty,
                new Binding() {Source = _osObject, Path = new PropertyPath("Name"), Mode = BindingMode.TwoWay});

            UpdateForeground();

            _osObject.PropertyChanged += (sender, args) =>
            {
                if (args.PropertyName.Equals("IsSelected"))
                {
                    if (IsSelected && !_osObject.IsSelected && !IsFocused)
                        IsSelected = false;
                    UpdateForeground();
                }
            };

            _osObject.ObjectChanged += (sender, args) =>
            {
                if (OpticalSceneObject.EventIDs.Disposed.Equals(args.ID))
                    Dispose();
            };
        }

        private void Dispose()
        {
            _tree.Items.Remove(this);
        }

        protected override void OnSelected(RoutedEventArgs e)
        {
            if (_tree != null)
                _tree.OSPresenter.ExclusiveSelect(_osObject);
            base.OnSelected(e);
        }

        internal void UpdateForeground()
        {
            if (_osObject == null) return;

            if (_osObject.IsSelected)
                Foreground = ResourceController.SelectedTextForeground;
            else
                Foreground = ResourceController.TextForeground;
        }
    }
}
