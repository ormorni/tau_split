# FastSpie 3

The FastSpie3 is similar to the FastSpie2 in the splitting method, but attempt not to do $O(R)$ work per step.

We will define two new terms:

* A reaction is **stable** when given the current bounds, the number of events is constant.
* A reaction is **active** as long as it is pushed down during the recursion.

We will deactivate a reaction when:

* It is stable.
* Either all reactions depending on its outputs are stable, or its event count is 0.

A question we now have to consider is when to reactivate a reaction. A reaction should be reactivated when its input product shift out of the range where is was confirmed to be stable. The most naïve method is to check for reactivation every time an input component of a reaction is updated. However, this is inefficient. A better way is to reactivate the reaction and split it as necessary only when there is a second-order reaction updating it. How do we determine that? An RData includes a second order reaction if when resampled, its event count changes. Only when that happens can it mess up the accounting of other reactions. When that happens, we split every reaction depending on it.

This is really pretty.

When does a reaction have to be active? When there is a dependent unstable reaction. I think it all works out.



## Deactivation and Reactivation

We need to take care when deactivating and reactivating reactions so as not to do excessive work.

There are two types of reactivation:

* Forward reactivation, where a reaction has second-order events that affect other reactions.
* Backward reactivations, where a formerly stable reaction was destabilized, and has to recompute its stability.

We will need data structures that allow both.

### Forward reactivation

The forward reactivation is simple, since it trivially occurs only when a reaction has a second-order event.

We will store a `Vec<Vec<ReactionData>>`, storing for every internal node the reaction data deactivated during that step, and for every leaf node the active reaction data.



### Backward reactivation

This one is more complex. We may need to reactivate all reaction data, but this has to be done step-by-step to make sure that the complexity analysis holds. We need for every component a vector storing the inactive reactions in each depth. We will then iterate over the reactions by iteration order and push them down step by step. We can store two data structure, one holding reactions by deactivation order, and one storing indices inside this vector showing which slice corresponds to which iteration node. Then we can easily push down inactive reactions.

Thus. we have a `Vec<Vec<usize>>` storing for each component the inactive reactions, and a `Vec<Vec<(u64, usize)>>`  storing the recursion indices of each block.



We will store a `Vec<Option<(usize, usize)>>`, pointing from every reaction to its position in the `Vec<Vec<ReactionData>>`, thus allowing quick access and removal of the inactive reaction data.



When iterating over the `by_component_stages` data structure, we note that while it is possible, I don't currently support random-access removal of elements from it, and the elements in it might not be up-to-date.

### Full interface

The data structure is:

```Rust
struct RecursionData {
    /// The current stages of the recursion.
    recursion_stages: Vec<u64>,
    /// Forward reactivation data structure
    inactive_by_stage: Vec<Vec<usize>>,
    /// Backward reactivation data structure
    inactive_by_component: Vec<Vec<Option<usize>>>,
    by_component_stages: Vec<Vec<(u64, usize)>>,
    /// The index mapping reactions to their position in the data structure.
    inactive_indices: Vec<Vec<Option<(usize, usize)>>>,
}

impl RecursionData {
    /// Adds a recursion step. 
    fn add_stage(&mut self, stage: u64);
    /// Removes a recursion step.
    /// Returns all inactive reactions remaining for that step.
    fn pop_stage(&mut self, stage: u64) -> Vec<ReactionData>;
    /// Completely reactivates a reaction.
    fn reactivate_reaction(&mut self, reaction: usize) -> ReactionData;
    /// Adds an inactive reaction to the data structure.
    fn add_inactive_reaction(&mut self, reaction_data: ReactionData);
    fn split_component(&mut self, component: usize);
}
```







## Complexity analysis





## Details

Now all that remains is to take care of the details.

At the beginning of the recursion, the active reactions are the last vector in `reaction_data`, the inactive are all non-leaf nodes, and the rest are just reactions that didn't yet take place. We need them in the stack so that we could reactivate reactions.

We then check reactivations. 

* We reactivate all reactions that were affected by a second-order step. These have to be updated, and maybe have a second-order reaction themselves.
* We reactivate all reactions required by a reaction that is no longer stable.





We then check if all reactions are stable. If they are, we apply all and end.

If not, we deactivate all reactions that can be deactivated and push them back to the stack. We then split the remaining reactions.

A minor optimization could be that an unstable reaction does not have to reactivate a zero reaction. I will do that later, after the prototype work.