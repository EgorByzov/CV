using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using Common;
using Common.TypeExtensions;
using ViewModelCore.Controllers;

namespace ViewModelCore.OpticalSceneEntities.SceneView
{
    internal class OSTreeView : TreeView
    {
        internal bool IsExclusiveSelect { get { return CurrentModif == ModifierKeys.None; } }
        internal bool IsExtraSelect { get { return CurrentModif == ModifierKeys.Shift; } }
        internal bool IsExpulsiveSelect { get { return CurrentModif == ModifierKeys.Control; } }

        internal OpenGLOpticalScene OSPresenter { get; private set; }

        private ModifierKeys CurrentModif { get; set; }

        protected override DependencyObject GetContainerForItemOverride()
        {
            return new OSTreeViewItem();
        }

        protected override bool IsItemItsOwnContainerOverride(object item)
        {
            return item is OSTreeViewItem;
        }

        protected override void OnPreviewMouseDown(MouseButtonEventArgs e)
        {
            CurrentModif = Keyboard.Modifiers;
            base.OnPreviewMouseDown(e);
        }

        protected override void OnPreviewKeyDown(KeyEventArgs e)
        {
            CurrentModif = Keyboard.Modifiers;
            base.OnPreviewKeyDown(e);
        }

        internal void InitCurrentProject(OpenGLOpticalScene osPresenter)
        {
            OSPresenter = osPresenter;
            Items.Clear();
            if (OSPresenter == null) return;

            var oso = OSPresenter.OS.SceneObjectTree;

            foreach (var elem in oso.Children)
                if (!elem.IsRegistrator)
                    Items.Add(new OSTreeViewItem(elem, this));

            OSPresenter.ObjectChanged += (sender, args) =>
            {
                var id = (OpenGLOpticalScene.EventIDs)args.ID;
                foreach (OSTreeViewItem item in Items)
                    item.OnOSEvent(id);
            };

            oso.ObjectChanged += (sender, args) =>
            {
                if (args.ID is OpticalSceneObject.ChildrenEventIDs)
                {
                    if (OpticalSceneObject.ChildrenEventIDs.Add.Equals(args.ID))
                    {
                        var obj = args.EventParams as OpticalSceneObject;
                        if (obj != null && !obj.IsRegistrator)
                            Items.Add(new OSTreeViewItem(obj, this));
                    }
                }
            };
        }

        internal void UpdateColors()
        {
            foreach (OSTreeViewItem item in Items)
                item.OnOSEvent(OpenGLOpticalScene.EventIDs.Root);
        }
    }

    public class OSTreeViewItem : TreeViewItem, INotifyPropertyChanged
    {
        #region Dependency

        private static readonly DependencyPropertyKey IsItemSelectedKey
        = DependencyProperty.RegisterReadOnly("IsItemSelected", typeof(bool), typeof(OSTreeViewItem), new FrameworkPropertyMetadata(false));
        public static readonly DependencyProperty IsItemSelectedProperty = IsItemSelectedKey.DependencyProperty;

        public bool IsItemSelected
        {
            get { return (bool)GetValue(IsItemSelectedProperty); }
            protected set { SetValue(IsItemSelectedKey, value); }
        }

        private static readonly DependencyPropertyKey IsRootExpandedKey
        = DependencyProperty.RegisterReadOnly("IsRootExpanded", typeof(bool), typeof(OSTreeViewItem), new FrameworkPropertyMetadata(false));
        public static readonly DependencyProperty IsRootExpandedProperty = IsItemSelectedKey.DependencyProperty;

        public bool IsRootExpanded
        {
            get { return (bool)GetValue(IsRootExpandedProperty); }
            protected set { SetValue(IsRootExpandedKey, value); }
        }

        private static readonly DependencyPropertyKey IsSelectableKey
        = DependencyProperty.RegisterReadOnly("IsSelectable", typeof(bool), typeof(OSTreeViewItem), new FrameworkPropertyMetadata(false));
        public static readonly DependencyProperty IsSelectableProperty = IsSelectableKey.DependencyProperty;

        public bool IsSelectable
        {
            get { return (bool)GetValue(IsSelectableProperty); }
            protected set { SetValue(IsSelectableKey, value); }
        }

        #endregion Dependency

        private OpticalSceneObject _osObject;
        private OSTreeView _tree;

        public event PropertyChangedEventHandler PropertyChanged;
        
        public String OSName { get { return _osObject.Name; } }

        // Flags
        public bool IsRegistrator { get { return _osObject.IsRegistrator; } }
        public bool IsSource { get { return _osObject.IsSource; } }
        public bool IsChecked
        {
            get { return _osObject.IsRoot || _osObject.IsRootAncestor; }
            set
            {
                if (value)
                    _tree.OSPresenter.SetCurrentRoot(_osObject);
                else if (_osObject.Parent != null)
                    _tree.OSPresenter.SetCurrentRoot(_osObject.Parent);
            }
        }
        public bool IsHidden { get { return _osObject.IsHidden; } set { _osObject.IsHidden = value; } }

        internal OSTreeViewItem() { }

        // TODO OS Object Property changed
        internal OSTreeViewItem(OpticalSceneObject osObj, OSTreeView tree)
        {
            _osObject = osObj;
            _tree = tree;

            UpdateDependencies();
            UpdateColor();


            _osObject.PropertyChanged += (sender, args) =>
            {
                switch (args.PropertyName)
                {
                    case "Name":
                        OnPropertyChanged("OSName");
                        break;
                    case "IsHidden":
                        OnPropertyChanged("IsHidden");
                        break;
                    default:
                        UpdateDependencies();
                        break;
                };
                UpdateColor();
            };

            foreach (var element in _osObject.Children)
                Items.Add(new OSTreeViewItem(element, tree));

            _osObject.ObjectChanged += (sender, args) =>
            {
                if (args.ID is OpticalSceneObject.ChildrenEventIDs)
                {
                    if (OpticalSceneObject.ChildrenEventIDs.Add.Equals(args.ID))
                    {
                        var obj = args.EventParams as OpticalSceneObject;
                        if (obj != null && !obj.IsRegistrator)
                            Items.Add(new OSTreeViewItem(obj, _tree));
                    }else if (OpticalSceneObject.ChildrenEventIDs.Remove.Equals(args.ID))
                    {
                        var obj = args.EventParams as OpticalSceneObject;

                    }
                }
                else if (args.ID is OpticalSceneObject.EventIDs && OpticalSceneObject.EventIDs.Disposed.Equals(args.ID) && Parent is ItemsControl)
                {
                    ((ItemsControl)Parent).Items.Remove(this);
                }
            };
        }

        protected override DependencyObject GetContainerForItemOverride()
        {
            return new OSTreeViewItem();
        }

        protected override bool IsItemItsOwnContainerOverride(object item)
        {
            return item is OSTreeViewItem;
        }

        internal void OnOSEvent(OpenGLOpticalScene.EventIDs eId)
        {
            if (eId == OpenGLOpticalScene.EventIDs.Root)
            {
                UpdateDependencies();
                OnPropertyChanged("IsChecked");
                UpdateColor();
            }

            foreach (OSTreeViewItem item in Items)
                item.OnOSEvent(eId);
        }

        private void UpdateDependencies()
        {
            IsRootExpanded = IsChecked;
            IsSelectable = _osObject.Parent != null && _osObject.Parent.IsRoot;
            IsItemSelected = _osObject.IsSelected;
        }

        protected override void OnMouseLeftButtonDown(MouseButtonEventArgs e)
        {
            if (!IsSelectable) return;

            if (_tree.IsExclusiveSelect)
                _tree.OSPresenter.ExclusiveSelect(_osObject);
            else if (_tree.IsExpulsiveSelect)
                _tree.OSPresenter.ExpulsiveSelect(_osObject);
            else if (_tree.IsExtraSelect)
                _tree.OSPresenter.ExtraSelect(_osObject);
        }

        protected override void OnMouseDoubleClick(MouseButtonEventArgs e)
        {
            if (!_osObject.IsRoot)
                _tree.OSPresenter.SetCurrentRoot(_osObject);
            else if (_osObject.Parent != null)
                _tree.OSPresenter.SetCurrentRoot(_osObject.Parent);
        }

        private void UpdateColor()
        {
            if (IsItemSelected)
                Foreground = ResourceController.SelectedTextForeground;
            else if (_osObject.IsRoot)
                Foreground = ResourceController.RootTextForeground;
            else if (!IsSelectable)
                Foreground = ResourceController.GhostedTextForeground;
            else
                Foreground = ResourceController.TextForeground;
        }

        protected virtual void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            var handler = PropertyChanged;
            if (handler != null) handler(this, new PropertyChangedEventArgs(propertyName));
        }
    }
}
