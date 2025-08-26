import json
import os
from typing import Union
from . import logger

def _to_internal_name(name: str) -> str:
    """
    For internal use. Convert a name to its internal (underscore-based) form.
    """
    return name.replace(" ", "_").replace("/", "_").replace("-", "_")

class TreeNode():
    """
    Class representing a single node in a cell type tree or hierarchy.

    Parameters
    ----------
    original_name
        The original label of the node (cell type), required and unique. This is what displays in the hierarchy.
    **kwargs
        Additional keyword arguments for describing the node, such as `cell_ontology_id`, `node_description`, `tissue_origin`, `markers`, `size`, `model`, and any other custom attributes.

    Attributes
    ----------
    original_name
        The original and display name of the node (cell type).
    internal_name
        A programmatic version of the original name with special characters replaced by underscores.
    cell_ontology_id
        Reference ID of the node (cell type) from the controlled Cell Ontology. Empty string if not provided.
    node_description
        Description of the node (cell type). Empty string if not provided.
    tissue_origin
        A list of tissue sources of the node (cell type). Empty list if not provided.
    markers
        A list of marker genes of the node (cell type). Empty list if not provided.
    size
        Number of cells contained in this node (cell type), typically populated/replaced during hierarchical model training. 0 if not provided.
    children
        A list of :class:`~celltypist.tree.TreeNode` instances representing child nodes (cell types). An empty list is initialized (thus a leaf node).
        For internal nodes, child nodes can be further added by :meth:`~celltypist.tree.TreeNode.add_children` after initialization, or directly from a json-like input (dict or json file).
        For leaf nodes, simply initialize them, or directly from a child-free json-like input (dict or json file).
    model
        The path to the CellTypist model used for classifying child cell types of the given internal node, typically populated/replaced during hierarchical model training. Empty string if not provided.
    """
    _STANDARD_FIELDS = ["cell_ontology_id", "node_description", "tissue_origin", "markers", "size", "children", "model"]

    def __init__(self, original_name: str, **kwargs):
        self.original_name = original_name
        if "internal_name" in kwargs:
            logger.warn("⚠️ The `internal_name` argument will be ignored because it is always derived from `original_name`")
        self.internal_name = _to_internal_name(original_name)
        self.cell_ontology_id = kwargs.get("cell_ontology_id", "")
        self.node_description = kwargs.get("node_description", "")
        self.tissue_origin = kwargs.get("tissue_origin", [])
        self.markers = kwargs.get("markers", [])
        self.size = kwargs.get("size", 0)
        if "children" in kwargs:
            logger.warn("⚠️ The `children` argument will be ignored. Use `add_children()` afterwards or read from a JSON structure instead")
        self.children = []
        self.model = kwargs.get("model", "")
        for key, val in kwargs.items():
            if key not in self._STANDARD_FIELDS and key != "internal_name":
                setattr(self, key, val)

    def add_children(self, *child_nodes) -> None:
        """
        Add/append child node(s) to the node's child list.

        Parameters
        ----------
        child_nodes
            One or more :class:`~celltypist.tree.TreeNode` instances.
            Can be passed individually or as a list/tuple/set.

        Returns
        ----------
        None
            Child node(s) are appended to the child list of the given node.
        """
        if not child_nodes:
            raise ValueError(
                    f"🛑 Please provide at least one child node to add to '{self.original_name}'")
        child_nodes = child_nodes[0] if isinstance(child_nodes[0], (list, tuple, set)) else child_nodes
        if not all(isinstance(child_node, TreeNode) for child_node in child_nodes):
            raise TypeError(
                    f"🛑 All child nodes to be added to '{self.original_name}' must be `TreeNode` instances")
        for child_node in child_nodes:
            if self.has_child(child_node.original_name):
                raise ValueError(
                        f"🛑 Duplicate child '{child_node.original_name}' should not be added to '{self.original_name}'")
        self.children.extend(child_nodes)

    def is_leaf(self) -> bool:
        """
        Check if this node is a leaf node.

        Returns
        ----------
        bool
            `True` if the node has no children, otherwise `False`.
        """
        return len(self.children) == 0

    def has_child(self, child_node) -> bool:
        """
        Check if a child node exists for the given node.

        Parameters
        ----------
        child_node
            A :class:`~celltypist.tree.TreeNode` instance or its `original_name` as a string.

        Returns
        ----------
        bool
            `True` if the child exists in the current node, otherwise `False`.
        """
        if isinstance(child_node, TreeNode):
            return child_node in self.children
        elif isinstance(child_node, str):
            return any(child.original_name == child_node for child in self.children)
        else:
            raise TypeError(
                    f"🛑 The query for '{self.original_name}' must be a `TreeNode` object or its original name (str)")

    def find_child(self, name: str):
        """
        Search for the child node with the given name.

        Parameters
        ----------
        name
            The original name of the child node to search for.

        Returns
        ----------
        :class:`~celltypist.tree.TreeNode` or None
            The matching child node if found; otherwise `None`.
        """
        return next((child for child in self.children if child.original_name == name), None)

    @property
    def child_names(self) -> list:
        """A list of original names of all immediate child nodes of the given node."""
        return [child.original_name for child in self.children]

    def remove_children(self, *child_nodes) -> None:
        """
        Remove one or more child nodes from the node's child list.

        Parameters
        ----------
        child_nodes
            One or more :class:`~celltypist.tree.TreeNode` instances, their `original_name`s as strings, or a mix of both.
            Can be passed individually or as a list/tuple/set.

        Returns
        ----------
        None
            Child node(s) are removed from the child list of the given node.
        """
        if not child_nodes:
            raise ValueError(
                    f"🛑 Please provide at least one child node to remove from '{self.original_name}'")
        child_nodes = child_nodes[0] if isinstance(child_nodes[0], (list, tuple, set)) else child_nodes
        if not all(isinstance(x, (TreeNode, str)) for x in child_nodes):
            raise TypeError(
                    f"🛑 All child nodes to be removed from '{self.original_name}' must be `TreeNode` instances or their original names (str)")
        for child_node in child_nodes:
            if not self.has_child(child_node):
                raise ValueError(
                        f"🛑 '{child_node.original_name if isinstance(child_node, TreeNode) else child_node}' is not a recognized child of '{self.original_name}'")
        to_remove = [x if isinstance(x, TreeNode) else next(child for child in self.children if child.original_name == x) for x in child_nodes]
        for child in to_remove:
            self.children.remove(child)

    def update(self, **kwargs) -> None:
        """
        Update attributes of the node with provided keyword arguments.

        Parameters
        ----------
        **kwargs
            Key-value pairs representing attribute names and their new values.

        Returns
        ----------
        None
            Attributes get updated.
        """
        for key, value in kwargs.items():
            if key in ('original_name', 'internal_name'):
                raise AttributeError(
                        f"🛑 '{key}' is read-only and cannot be updated for '{self.original_name}'")
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise AttributeError(
                        f"🛑 '{key}' is not a valid attribute for updating '{self.original_name}'")
        self.validate(check_type = True)

    def validate(self, check_type: bool = True) -> None:
        """
        Sanity check the node.

        Parameters
        ----------
        check_type
            Whether to validate types of the node's attributes.
            (Default: `True`)

        Returns
        ----------
        None
            The given node is checked.
        """
        if check_type:
            for attr in ("original_name", "cell_ontology_id", "node_description", "model"):
                if not isinstance(getattr(self, attr), str):
                    raise TypeError(
                            f"🛑 `{attr}` of the node '{self.original_name}' must be a string")
            for attr in ("tissue_origin", "markers"):
                val = getattr(self, attr)
                if not isinstance(val, list) or not all(isinstance(x, str) for x in val):
                    raise TypeError(
                            f"🛑 `{attr}` of the node '{self.original_name}' must be a list of strings")
            if not isinstance(self.size, int):
                raise TypeError(
                        f"🛑 `size` of the node '{self.original_name}' must be an integer")
            if not isinstance(self.children, list) or not all(isinstance(c, TreeNode) for c in self.children):
                raise TypeError(
                        f"🛑 `children` of the node '{self.original_name}' must be a list of `TreeNode` instances")
        if self.is_leaf():
            if self.model:
                raise ValueError(
                        f"🛑 Leaf node '{self.original_name}' should not have a model")
            return
        if len(set(self.child_names)) != len(self.child_names):
            raise ValueError(
                    f"🛑 Duplicate child names found under the node '{self.original_name}'")
        if self.original_name in self.child_names:
            raise ValueError(
                    f"🛑 Node '{self.original_name}' cannot have a child with the same name")
        for child in self.children:
            child.validate(check_type = check_type)

    def to_dict(self) -> dict:
        """
        Convert the node and its children into a dictionary representation.

        Returns
        ----------
        dict
            A dictionary with the node's non-empty attributes.
            Child nodes (if any) are recursively included under the `children` key.
        """
        data = {}
        for key, value in self.__dict__.items():
            if key == "internal_name" or (key == "size" and value == 0):
                continue
            if key == "children":
                value = [child.to_dict() for child in self.children]
            if value not in (None, '', []):
                data[key] = value
        return data

    def write(self, file_path: str) -> None:
        """Write the node to a JSON file."""
        file_path = os.path.splitext(file_path)[0] + '.json'
        with open(file_path, "w") as f:
            json.dump(self.to_dict(), f, indent = 2)

    @classmethod
    def from_dict(cls, node_dict: dict):
        """Build a :class:`~celltypist.tree.TreeNode` instance from a dictionary."""
        if "original_name" not in node_dict:
            raise KeyError(
                    f"🛑 Each node must have an 'original_name'")
        original_name = node_dict["original_name"]
        node_data = {k: v for k, v in node_dict.items() if k not in ("original_name", "internal_name", "children")}
        node = cls(original_name, **node_data)
        node.children = [cls.from_dict(child) for child in node_dict.get("children", [])]
        node.validate(check_type = True)
        return node

    @classmethod
    def from_json(cls, json_data: Union[dict, str]):
        """
        Build a :class:`~celltypist.tree.TreeNode` instance from a dictionary or JSON file.

        Parameters
        ----------
        json_data
            Either a dictionary representing a serialized :class:`~celltypist.tree.TreeNode`, or a file path to a JSON file containing such a dictionary.

        Returns
        ----------
        :class:`~celltypist.tree.TreeNode`
            A reconstructed :class:`~celltypist.tree.TreeNode` instance with full hierarchy.
        """
        if isinstance(json_data, str):
            with open(json_data, "r") as f:
                json_data = json.load(f)
        return cls.from_dict(json_data)

    def __repr__(self):
        """String representation of the :class:`~celltypist.tree.TreeNode` object."""
        base = f"A tree node '{self.original_name}' with {len(self.children)} child nodes"
        for x in self._STANDARD_FIELDS:
            val = getattr(self, x)
            if val not in ('', [], 0):
                if x == 'children':
                    val = self.child_names
                if x in ['tissue_origin', 'markers', 'children']:
                    if len(val) == 1:
                        val = val[0]
                    elif len(val) == 2:
                        val = f"{val[0]}, {val[1]}"
                    elif len(val) == 3:
                        val = f"{val[0]}, {val[1]}, {val[2]}"
                    else:
                        val = f"{val[0]}, {val[1]}, ..., {val[-1]}"
                base += f"\n    {x}: {val}"
        return base
