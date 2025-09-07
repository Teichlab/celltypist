import json
import os
import copy
from typing import Union, Optional, Iterator
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
        For internal nodes, child nodes can be further added by :meth:`~celltypist.tree.TreeNode.add_children` after initialization, or directly from a JSON-like input (dict or JSON file).
        For leaf nodes, simply initialize them, or directly from a child-free JSON-like input (dict or JSON file).
    model
        The path to the CellTypist model used for classifying child cell types of the given internal node, typically populated/replaced during hierarchical model training. Empty string if not provided.
    child_names
        A list of child original names for the node. Empty list for a leaf node.
    depth
        The depth of the node. For example, a node with three leaf children has a depth of 2.
    n_leaves
        The number of leaf nodes contained.
    n_nodes
        The number of total nodes contained.
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

    def add_children(self, *child_nodes) -> list:
        """
        Add/append child node(s) to the node's child list.

        Parameters
        ----------
        child_nodes
            One or more :class:`~celltypist.tree.TreeNode` instances.
            Can be passed individually or as a list/tuple/set.

        Returns
        ----------
        list
            A list of child node(s) that were appended to the child list of the given node.
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
        return list(child_nodes)

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
            return any(child.original_name == child_node.original_name for child in self.children)
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

    @property
    def depth(self) -> int:
        """The depth of the node."""
        if self.is_leaf():
            return 1
        return 1 + max(child.depth for child in self.children)

    @property
    def n_leaves(self) -> int:
        """Number of leaves recursively contained."""
        if self.is_leaf():
            return 1
        return sum(child.n_leaves for child in self.children)

    @property
    def n_nodes(self) -> int:
        """Number of total nodes recursively contained."""
        return 1 + sum(child.n_nodes for child in self.children)

    def cell_types(self, leaf_only: bool = True) -> list:
        """
        Return cell types recursively contained in this node.

        Parameters
        ----------
        leaf_only
            Whether to return only leaf (finest-grained) or all (internal + leaf) cell types.
            (Default: `True`)

        Returns
        ----------
        list
            List of cell type names from this node.
        """
        names = []
        if not leaf_only or self.is_leaf():
            names.append(self.original_name)
        for child in self.children:
            names.extend(child.cell_types(leaf_only = leaf_only))
        return names

    def remove_children(self, *child_nodes) -> list:
        """
        Remove one or more child nodes from the node's child list.

        Parameters
        ----------
        child_nodes
            One or more :class:`~celltypist.tree.TreeNode` instances, their `original_name`s as strings, or a mix of both.
            Can be passed individually or as a list/tuple/set.

        Returns
        ----------
        list
            A list of child node(s) that were removed from the child list of the given node.
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
        return to_remove

    def reorder_children(self, new_order: Union[list, tuple]) -> list:
        """
        Reorder the direct children of the node to match `new_order`.

        Parameters
        ----------
        new_order
            List of child names specifying the desired order.

        Returns
        ----------
        list
            The reordered list of child nodes.
        """
        current = self.child_names
        if not current:
            return current
        if set(current) != set(new_order) or len(current) != len(new_order):
            raise ValueError(
                    f"🛑 `new_order` must be a permutation of the current child names for '{self.original_name}'")
        name_to_node = {c.original_name: c for c in self.children}
        self.children = [name_to_node[n] for n in new_order]
        return self.children

    def sort_children(self, recursive: bool = True, descending: bool = True) -> list:
        """
        Sort this node's children by the number of leaves they contain.

        Parameters
        ----------
        recursive
            Sort by the total number of leaves recursively contained (`recursive = True`), or by the number of direct children (`recursive = False`).
            (Default: `True`)
        descending
            Whether to sort in descending order.
            (Default: `True`)

        Returns
        ----------
        list
            The sorted list of child nodes.
        """
        if self.is_leaf():
            return []
        key_fn = (lambda c: c.n_leaves) if recursive else (lambda c: len(c.children))
        self.children.sort(key = key_fn, reverse = descending)
        return self.children

    def update(self, validate: bool = False, **kwargs) -> None:
        """
        Update attributes of the node with provided keyword arguments.

        Parameters
        ----------
        validate
            Whether to validate `kwargs` (keyword arguments) before updating.
            (Default: `False`)
        **kwargs
            Key-value pairs representing attribute names and their new values.

        Returns
        ----------
        None
            Attributes get updated.
        """
        for key in kwargs:
            if key == 'internal_name':
                raise AttributeError(
                        f"🛑 '{key}' cannot be updated for '{self.original_name}' as it will always be derived")
            if not hasattr(self, key):
                raise AttributeError(
                        f"🛑 '{key}' is not a valid attribute for updating '{self.original_name}'")
        if validate:
            proposed = copy.copy(self)
            for key, value in kwargs.items():
                setattr(proposed, key, value)
            proposed.validate(check_type = True, check_unique = True)
        for key, value in kwargs.items():
            setattr(self, key, value)
            if key == 'original_name':
                self.internal_name = _to_internal_name(value)

    def validate(self, check_type: bool = True, check_unique: bool = True) -> None:
        """
        Sanity check the node.

        Parameters
        ----------
        check_type
            Whether to validate types of the node's attributes.
            (Default: `True`)
        check_unique
            Whether to check uniqueness of cell type names recursively contained in the node.
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
        if self.is_leaf() and self.model:
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
            child.validate(check_type = check_type, check_unique = False)
        if check_unique:
            names = self.cell_types(leaf_only = False)
            duplicates = [n for i, n in enumerate(names) if n in names[:i]]
            if duplicates:
                raise ValueError(
                        f"🛑 Duplicate node names found: {sorted(set(duplicates))}")

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
        if not isinstance(node_dict, dict):
            raise TypeError(
                    f"🛑 Please provide a top-level JSON object (dict) for `from_dict`")
        if "original_name" not in node_dict:
            raise KeyError(
                    f"🛑 Each node must have an 'original_name'")
        original_name = node_dict["original_name"]
        node_data = {k: v for k, v in node_dict.items() if k not in ("original_name", "internal_name", "children")}
        node = cls(original_name, **node_data)
        node.children = [cls.from_dict(child) for child in node_dict.get("children", [])]
        node.validate(check_type = True, check_unique = True)
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

    def copy(self):
        """
        Create a deep copy of the node.

        Returns
        ----------
        :class:`~celltypist.tree.TreeNode`
            A deep-copied :class:`~celltypist.tree.TreeNode` instance.
        """
        return TreeNode.from_dict(self.to_dict())

    def __contains__(self, name: str):
        """Return True if a node with the given name exists."""
        if not isinstance(name, str):
            raise TypeError(
                    f"🛑 Membership check requires a string name.")
        return name in self.cell_types(leaf_only = False)

    def __repr__(self):
        """String representation of the :class:`~celltypist.tree.TreeNode` object."""
        n = len(self.children)
        if n == 0:
            base = f"A leaf node '{self.original_name}'"
        else:
            base = f"An internal node '{self.original_name}' with {n} child {'node' if n == 1 else 'nodes'}"
        for x in self._STANDARD_FIELDS:
            val = getattr(self, x)
            if val not in ('', [], 0):
                if x == 'children':
                    val = self.child_names
                if x in ('tissue_origin', 'markers', 'children'):
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

class Tree():
    """
    Tree-level wrapper around a :class:`~celltypist.tree.TreeNode` root.

    Parameters
    ----------
    handle
        Machine-friendly unique identifier/name of this cell type tree/hierarchy.
    root
        A :class:`~celltypist.tree.TreeNode` root of the tree.
    **kwargs
        Additional keyword arguments for describing the tree.

    Attributes
    ----------
    handle
        Unique programmatic identifier/name of this tree.
    root
        The root node (a :class:`~celltypist.tree.TreeNode` instance) of the tree.
    depth
        The depth of the tree.
    n_leaves
        The number of leaf nodes contained.
    n_nodes
        The number of total nodes contained.
    """
    def __init__(self, handle: str, root: TreeNode, **kwargs):
        if not isinstance(handle, str):
            raise TypeError(
                    f"🛑 `handle` must be a string")
        if handle.strip() == '':
            raise ValueError(
                    f"🛑 `handle` must be a non-empty string")
        if not isinstance(root, TreeNode):
            raise TypeError(
                    f"🛑 `root` must be a `TreeNode` instance")
        self.handle = handle
        self.root = root
        for key, value in kwargs.items():
            setattr(self, key, value)

    #def __getattr__(self, name):
    #    """
    #    Delegate unknown attribute/method lookups to the root node after standard MRO.
    #    """
    #    return getattr(self.root, name)

    def validate(self, check_type: bool = True, check_unique: bool = True) -> None:
        if not isinstance(self.handle, str):
            raise TypeError(
                    f"🛑 `handle` must be a string")
        if self.handle.strip() == '':
            raise ValueError(
                    f"🛑 `handle` must be a non-empty string")
        if not isinstance(self.root, TreeNode):
            raise TypeError(
                    f"🛑 `root` must be a `TreeNode` instance")
        self.root.validate(check_type = check_type, check_unique = check_unique)

    @property
    def depth(self) -> int:
        return self.root.depth

    @property
    def n_leaves(self) -> int:
        return self.root.n_leaves

    @property
    def n_nodes(self) -> int:
        return self.root.n_nodes

    def cell_types(self, leaf_only: bool = True) -> list:
        return self.root.cell_types(leaf_only = leaf_only)

    def to_dict(self) -> dict:
        """
        Convert the tree into a dictionary representation.

        Returns
        ----------
        dict
            Dictionary containing `handle` (unique identifier), `root` (serialized root node), and any custom attributes of the tree.
        """
        data = {"handle": self.handle}
        for key, value in self.__dict__.items():
            if key not in ("handle", "root"):
                data[key] = value
        data["root"] = self.root.to_dict()
        return data

    def write(self, file_path: str) -> None:
        """Write the tree to a JSON file."""
        file_path = os.path.splitext(file_path)[0] + '.json'
        with open(file_path, "w") as f:
            json.dump(self.to_dict(), f, indent = 2)

    @classmethod
    def from_dict(cls, tree_dict: dict):
        """Build a :class:`~celltypist.tree.Tree` instance from a dictionary."""
        if not isinstance(tree_dict, dict):
            raise TypeError(
                    f"🛑 Please provide a top-level JSON object (dict) for `from_dict`")
        if "handle" not in tree_dict:
            raise KeyError(
                    f"🛑 Each tree must have a 'handle'")
        if "root" not in tree_dict:
            raise KeyError(
                    f"🛑 Each tree must have a 'root'")
        extras = {k: v for k, v in tree_dict.items() if k not in ("handle", "root")}
        return cls(handle = tree_dict["handle"], root = TreeNode.from_dict(tree_dict["root"]), **extras)

    @classmethod
    def from_json(cls, json_data: Union[dict, str]):
        if isinstance(json_data, str):
            with open(json_data, "r") as f:
                json_data = json.load(f)
        return cls.from_dict(json_data)

    def copy(self):
        return Tree.from_dict(self.to_dict())

    def __contains__(self, name: str):
        """Return True if a node with the given name exists."""
        return name in self.root

    def __repr__(self):
        """String representation of the :class:`~celltypist.tree.Tree` object."""
        base = f"A cell type tree with {self.n_nodes} total nodes and {self.n_leaves} leaves"
        base += f"\n    handle: {self.handle}"
        for key, value in self.__dict__.items():
            if key not in ("handle", "root"):
                base += f"\n    {key}: {value}"
        base += f"\n    root: a node '{self.root.original_name}' with depth {self.root.depth}"
        return base

    @staticmethod
    def _traverse(node: TreeNode, parent: Optional[TreeNode] = None) -> Iterator:
        """
        For internal use. Yield (node, parent) pairs in a DFS traversal.
        """
        #if isinstance(parent, TreeNode):
        #    if not parent.has_child(node):
        #        raise ValueError(
        #                f"🛑 Please provide matched `node` and `parent`")
        yield node, parent
        for child in node.children:
            yield from Tree._traverse(child, parent = node)

    def iter_nodes(self, leaf_only: bool = False) -> Iterator:
        """A generator iterating through the tree."""
        for node, _ in Tree._traverse(self.root):
            if leaf_only:
                if node.is_leaf():
                    yield node
            else:
                yield node

    def find_node(self, name: str) -> TreeNode:
        """
        Find a node in the tree by its name.

        Parameters
        ----------
        name
            The node name to search.

        Returns
        ----------
        :class:`~celltypist.tree.TreeNode`
            A :class:`~celltypist.tree.TreeNode` instance with the given name.
        """
        if not isinstance(name, str):
            raise TypeError(
                    f"🛑 `name` must be a string")
        for node, _ in Tree._traverse(self.root):
            if node.original_name == name:
                return node
        raise ValueError(
                f"🛑 No node named '{name}' exists in this tree")

    def find_parent(self, name: str) -> TreeNode:
        """
        Find the parent of a node in the tree by the node's name.

        Parameters
        ----------
        name
            The name of the node whose parent is queried.

        Returns
        ----------
        :class:`~celltypist.tree.TreeNode`
            A :class:`~celltypist.tree.TreeNode` instance representing the parent of the specified node.
        """
        if not isinstance(name, str):
            raise TypeError(
                    f"🛑 `name` must be a string")
        for node, parent in Tree._traverse(self.root):
            if node.original_name == name:
                if parent is None:
                    raise ValueError(
                            f"🛑 Node '{name}' is the root and has no parent")
                return parent
        raise ValueError(
                f"🛑 No node named '{name}' exists in this tree")

    def find_siblings(self, name: str, return_names: bool = True) -> list:
        """
        Return siblings of a node (nodes with the same parent), excluding the node itself.

        Parameters
        ----------
        name
            The name of the node whose siblings are queried.
        return_names
            Whether to return a list of sibling names (rather than a list of :class:`~celltypist.tree.TreeNode` objects).
            (Default: `True`)

        Returns
        ----------
        list
            A list of sibling names (if `return_names = True`) or :class:`~celltypist.tree.TreeNode` objects (if `return_names = False`).
        """
        if not isinstance(name, str):
            raise TypeError(
                    f"🛑 `name` must be a string")
        if self.root.original_name == name:
            raise ValueError(
                    f"🛑 Node '{name}' is the root and has no siblings")
        parent = self.find_parent(name)
        nodes = [c for c in parent.children if c.original_name != name]
        return [n.original_name for n in nodes] if return_names else nodes

    def find_children(self, name: str, return_names: bool = True) -> list:
        """
        Find the children of a node in the tree by the node's name.

        Parameters
        ----------
        name
            The name of the node whose children are queried.
        return_names
            Whether to return a list of child names (rather than a list of :class:`~celltypist.tree.TreeNode` objects).
            (Default: `True`)

        Returns
        ----------
        list
            A list of child names (if `return_names = True`) or :class:`~celltypist.tree.TreeNode` objects (if `return_names = False`).
        """
        node = self.find_node(name)
        #node.children is not protected; any modification is MIP.
        return node.child_names if return_names else node.children

    def add_children(self, parent: str, *child_nodes) -> list:
        """
        Add one or more children to a node in the tree by the node's name.

        Parameters
        ----------
        parent
            The name of the parent node to which children are added.
        *child_nodes
            One or more :class:`~celltypist.tree.TreeNode` instances.
            Can be passed individually or as a list/tuple/set.

        Returns
        ----------
        list
            A list of child node(s) that were appended to the child list of the given node.
        """
        if not isinstance(parent, str):
            raise TypeError(
                    f"🛑 `parent` must be a string")
        node = self.find_node(parent)
        return node.add_children(*child_nodes)

    def remove_children(self, parent: str, *child_nodes) -> list:
        """
        Remove one or more children from a node in the tree by the node's name.

        Parameters
        ----------
        parent
            The name of the parent node whose children are to be removed.
        child_nodes
            One or more :class:`~celltypist.tree.TreeNode` instances, their `original_name`s as strings, or a mix of both.
            Can be passed individually or as a list/tuple/set.

        Returns
        ----------
        list
            A list of child node(s) that were removed from the child list of the given node.
        """
        if not isinstance(parent, str):
            raise TypeError(
                    f"🛑 `parent` must be a string")
        node = self.find_node(parent)
        return node.remove_children(*child_nodes)

    def reorder_children(self, parent: str, new_order: Union[list, tuple]) -> list:
        """
        Reorder the direct children of a parent node to match `new_order`.

        Parameters
        ----------
        parent
            The name of the parent node whose children will be reordered.
        new_order
            List of child names specifying the desired order.

        Returns
        ----------
        list
            The reordered list of child nodes.
        """
        if not isinstance(parent, str):
            raise TypeError(
                    f"🛑 `parent` must be a string")
        node = self.find_node(parent)
        return node.reorder_children(new_order)

    def sort_children(self, parent: str, recursive: bool = True, descending: bool = True) -> list:
        """
        Sort the direct children of a parent node by the number of leaves they contain.

        Parameters
        ----------
        parent
            The name of the parent node whose children will be sorted.
        recursive
            Sort by the total number of leaves recursively contained (`recursive = True`), or by the number of direct children (`recursive = False`).
            (Default: `True`)
        descending
            Whether to sort in descending order.
            (Default: `True`)

        Returns
        ----------
        list
            The sorted list of child nodes.
        """
        if not isinstance(parent, str):
            raise TypeError(
                    f"🛑 `parent` must be a string")
        node = self.find_node(parent)
        return node.sort_children(recursive = recursive, descending = descending)

    @staticmethod
    def _compute_n_leaves(node: TreeNode, counts: dict) -> int:
        """
        For internal use. Populate `counts` with entries for each node and returns the leaf count of the current node.
        """
        if node.is_leaf():
            counts[node.original_name] = 1
        else:
            counts[node.original_name] = sum(Tree._compute_n_leaves(c, counts) for c in node.children)
        return counts[node.original_name]

    @property
    def n_leaves_by_node(self) -> dict:
        """Return a mapping. Keys are node names, values are total leaf counts recursively contained under that node."""
        counts = {}
        Tree._compute_n_leaves(self.root, counts)
        return counts

    def remove_node(self, name: str) -> TreeNode:
        """
        Remove a node from the tree by its name.

        Parameters
        ----------
        name
            The name of the node to remove.

        Returns
        ----------
        :class:`~celltypist.tree.TreeNode`
            The removed :class:`~celltypist.tree.TreeNode` instance.
        """
        if not isinstance(name, str):
            raise TypeError(
                    f"🛑 `name` must be a string")
        if self.root.original_name == name:
            raise ValueError(
                    f"🛑 Cannot remove the root node from the tree")
        parent = self.find_parent(name)
        return parent.remove_children(name)[0]

    def add_node(self, parent: str, node: TreeNode) -> TreeNode:
        """
        Add a node under a given parent node in the tree.

        Parameters
        ----------
        parent
            The name of the parent node to which the new node will be added.
        node
            The :class:`~celltypist.tree.TreeNode` instance to add.

        Returns
        ----------
        :class:`~celltypist.tree.TreeNode`
            The :class:`~celltypist.tree.TreeNode` instance that was added.
        """
        return self.add_children(parent, node)[0]

    def move_node(self, name: str, to: str, validate: bool = True) -> TreeNode:
        """
        Move a node (and its subtree) from its current parent to a new parent.

        Parameters
        ----------
        name
            The name of the node to move.
        to
            The name of the new parent node.
        validate
            Whether to validate `name` and `to`.
            (Default: `True`)

        Returns
        ----------
        :class:`~celltypist.tree.TreeNode`
            The moved :class:`~celltypist.tree.TreeNode` instance.
        """
        if not isinstance(name, str) or not isinstance(to, str):
            raise TypeError(
                    f"🛑 Both `name` and `to` must be strings")
        if self.root.original_name == name:
            raise ValueError(
                    f"🛑 Cannot move the root node")
        if validate:
            cell_types = self.cell_types(leaf_only = False)
            if name not in cell_types:
                raise ValueError(
                        f"🛑 No node named '{name}' exists in this tree")
            if to not in cell_types:
                raise ValueError(
                        f"🛑 No node named '{to}' exists in this tree")
        removed = self.remove_node(name)
        return self.add_node(to, removed)

    def update(self, name: str, validate: bool = False, **kwargs) -> None:
        """
        Update attributes of a node in the tree by the node's name.

        Parameters
        ----------
        name
            The name of the node to update.
        validate
            Whether to validate the updated node and its descendants after applying the update.
            (Default: `False`)
        **kwargs
            Key-value pairs of attributes to update on the node.

        Returns
        ----------
        None
            Attributes of the given node get updated.
        """
        node = self.find_node(name)
        node.update(validate = validate, **kwargs)

    def replace_node(self, name: str, by: TreeNode) -> TreeNode:
        """
        Replace a node in the tree (and its subtree) with a new node.

        Parameters
        ----------
        name
            The name of the node to replace.
        by
            The new :class:`~celltypist.tree.TreeNode` to replace.

        Returns
        ----------
        :class:`~celltypist.tree.TreeNode`
            The old :class:`~celltypist.tree.TreeNode` that was replaced.
        """
        if not isinstance(name, str):
            raise TypeError(
                    f"🛑 `name` must be a string")
        if self.root.original_name == name:
            raise ValueError(
                    f"🛑 Cannot replace the root node")
        if name not in self:
            raise ValueError(
                    f"🛑 '{name}' does not exist in the tree")
        if not isinstance(by, TreeNode):
            raise TypeError(
                    f"🛑 `by` must be a `TreeNode`")
        parent = self.find_parent(name)
        parent.add_children(by)
        return parent.remove_children(name)[0]

    def extract_subtree(self, name: str, handle: Optional[str] = None):
        """
        Extract a subtree rooted at a given node, returning it as a new tree.

        Parameters
        ----------
        name
            The name of the node to use as the root of the new subtree.
        handle
            Handle for the new tree. Defaults to `f"subtree_of_{name}"`.

        Returns
        ----------
        :class:`~celltypist.tree.Tree`
            A new :class:`~celltypist.tree.Tree` instance rooted at the specified node.
        """
        proposed = self.copy()
        proposed.root = proposed.find_node(name)
        proposed.handle = handle or f"subtree_of_{name}"
        return proposed

    def extract_path(self, name: str) -> list:
        """
        Extract the path from the root to a given node.

        Parameters
        ----------
        name
            The name of the target node.

        Returns
        ----------
        list
            A list of nodes along the path from root to the target.
        """
        if not isinstance(name, str):
            raise TypeError(
                    f"🛑 `name` must be a string")
        if name not in self:
            raise ValueError(
                    f"🛑 '{name}' does not exist in the tree")
        parent_map = {}
        for node, parent in Tree._traverse(self.root):
            parent_map[node] = parent
            if node.original_name == name:
                break
        path = []
        while node is not None:
            path.append(node)
            node = parent_map[node]
        path.reverse()
        print(" -> ".join(n.original_name for n in path))
        return path

    def lowest_common_ancestor(self, name1: str, name2: str) -> TreeNode:
        """
        Find the lowest common ancestor of two nodes in the tree.

        Parameters
        ----------
        name1
            The name of the first node.
        name2
            The name of the second node.

        Returns
        ----------
        :class:`~celltypist.tree.TreeNode`
            The lowest common ancestor node of the two nodes.
        """
        path1 = self.extract_path(name1)
        path2 = self.extract_path(name2)
        for n1, n2 in zip(path1, path2):
            if n1 == n2:
                lca = n1
            else:
                break
        return lca

Tree.validate.__doc__ = TreeNode.validate.__doc__.replace("node", "tree")
Tree.depth.__doc__ = TreeNode.depth.__doc__.replace("node", "tree")
Tree.n_leaves.__doc__ = TreeNode.n_leaves.__doc__
Tree.n_nodes.__doc__ = TreeNode.n_nodes.__doc__
Tree.cell_types.__doc__ = TreeNode.cell_types.__doc__.replace("node", "tree")
Tree.from_json.__func__.__doc__ = TreeNode.from_json.__func__.__doc__.replace("TreeNode", "Tree")
Tree.copy.__doc__ = TreeNode.copy.__doc__.replace("node", "tree").replace("TreeNode", "Tree")
