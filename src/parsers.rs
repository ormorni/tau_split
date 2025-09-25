use crate::reaction::Reaction;
use derive_new::new;
use itertools::Itertools;
use nom::{
    branch::alt,
    bytes::complete::{tag, take_while1},
    character::complete::{digit0, digit1, multispace0},
    combinator::map_res,
    multi::separated_list0,
    number::complete::double,
    sequence::delimited,
    AsChar, IResult, Parser,
};
use num_traits::Zero;
use rustc_hash::FxHashMap;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};
use tinyvec::ArrayVec;

/// The result of parsing a line defining a new reaction.
#[derive(Clone, Debug)]
struct NamedReaction {
    inputs: Vec<(String, u64)>,
    outputs: Vec<(String, u64)>,
    rate: f64,
}

/// The result of parsing a line setting the initial value of a reactant.
#[derive(Clone, Debug, new)]
struct Reactant {
    name: String,
    amount: u64,
}

#[derive(Default)]
pub struct ParseState {
    initial_states: FxHashMap<String, u64>,
    reactions: Vec<NamedReaction>,
}

/// An enum storing the result of parsing a line.
enum Line {
    Reactant(Reactant),
    Reaction(NamedReaction),
}

/// A parser for a nonnegative decimal number.
pub fn decimal(data: &str) -> IResult<&str, u64> {
    map_res(digit1, |s: &str| s.parse::<u64>()).parse(data)
}

/// Parses a line of the form:
/// ```ignore
/// A = 5
/// ```
/// that sets the initial value of the reactant `A` to 5.
fn parse_reactant(data: &str) -> IResult<&str, Line> {
    let (rem, (name, _, _, _, amount)) = (
        take_while1(AsChar::is_alphanum),
        multispace0,
        tag("="),
        multispace0,
        decimal,
    )
        .parse(data)?;

    Ok((rem, Line::Reactant(Reactant::new(name.to_owned(), amount))))
}

/// Parses a term of the form `2A`.
fn parse_reaction_item(data: &str) -> IResult<&str, (String, u64)> {
    let (rem, num): (&str, u64) = map_res(digit0, |s: &str| {
        if s.is_empty() {
            Ok(1)
        } else {
            s.parse::<u64>()
        }
    })
    .parse(data)?;
    let (rem, name) = take_while1(AsChar::is_alphanum).parse(rem)?;

    Ok((rem, (name.to_owned(), num)))
}

/// Parses one-half of a reaction:
/// ```ignore
/// 2A + B
/// ```
fn parse_reaction_half(data: &str) -> IResult<&str, Vec<(String, u64)>> {
    separated_list0(
        delimited(multispace0, tag("+"), multispace0),
        parse_reaction_item,
    )
    .parse(data)
}

/// A parser for a full reaction, of the form:
/// ```ignore
/// 2A + B -> 3C, 3.5e-9
/// ```
fn parse_reaction(reaction: &str) -> IResult<&str, Line> {
    let (rem, (left_half, _, right_half, _, rate)) = (
        parse_reaction_half,
        delimited(multispace0, tag("->"), multispace0),
        parse_reaction_half,
        delimited(multispace0, tag(","), multispace0),
        double,
    )
        .parse(reaction)?;

    let res = NamedReaction {
        inputs: left_half,
        outputs: right_half,
        rate,
    };

    Ok((rem, Line::Reaction(res)))
}

fn parse_line(line: &str) -> IResult<&str, Line> {
    alt((parse_reactant, parse_reaction)).parse(line)
}

fn named_to_reaction(
    named_reaction: NamedReaction,
    reactant_names: &FxHashMap<String, usize>,
) -> Reaction {
    let mut inputs = ArrayVec::new();
    for (comp, count) in &named_reaction.inputs {
        let comp = *reactant_names.get(comp).unwrap_or_else(||panic!("Failed to parse the reaction: \"{named_reaction:?}\": The reactant \"{comp:?}\" is undefined!"));
        if inputs
            .last()
            .is_some_and(|&(last_comp, _)| last_comp == comp)
        {
            inputs.last_mut().unwrap().1 += *count;
        } else {
            inputs.push((comp, *count));
        }
    }

    let outputs = named_reaction
        .outputs
        .into_iter()
        .map(|(s, count)| {
            let comp = *reactant_names.get(&s).expect("Failed to parse the reaction: \"{named_reaction:?}\": The reactant \"{comp:?}\" is undefined!");
            (comp, count as i64)
        })
        .collect_vec();

    // Computing an iterator over the differences, and merging it to a single stoichiometry vector.
    let in_diff = inputs.iter().map(|(idx, count)| (*idx, -(*count as i64)));
    let all_diff = in_diff.chain(outputs.into_iter()).sorted();
    let mut stoichiometry: ArrayVec<[(usize, i64); 4]> = ArrayVec::new();
    for (idx, diff) in all_diff {
        if stoichiometry.is_empty() || stoichiometry.last().unwrap().0 < idx {
            stoichiometry.push((idx, diff));
        } else {
            stoichiometry.last_mut().unwrap().1 += diff;
        }
    }
    let stoichiometry = stoichiometry
        .iter()
        .filter(|(_, diff)| !diff.is_zero())
        .copied()
        .collect();

    Reaction::new(inputs, stoichiometry, named_reaction.rate)
}

impl ParseState {
    /// Parses a data file.
    /// The data file contains lines, each of which is either a definition of the initial state of a reactant:
    /// ```ignore
    /// A = 5
    /// B = 7
    /// ```
    /// or a reaction
    pub fn parse_data_file(&mut self, reactions_path: &Path) -> &mut Self {
        BufReader::new(
            File::open(reactions_path)
                .unwrap_or_else(|err| panic!("Failed to open {reactions_path:?}: {err:?}!")),
        )
        .lines()
        .filter_map(|line| line.ok())
        .filter(|line| !line.trim().starts_with("#"))
        .map(|line| {
            parse_line(&line)
                .unwrap_or_else(|err| panic!("Failed to parse the line {line} with error {err:?}"))
                .1
        })
        .for_each(|line| match line {
            Line::Reactant(reactant) => {
                self.initial_states.insert(reactant.name, reactant.amount);
            }
            Line::Reaction(named_reaction) => {
                self.reactions.push(named_reaction);
            }
        });
        self
    }

    /// Gets the reaction network.
    /// The network has three components:
    /// * The initial state.
    /// * The reactions.
    /// * The name of each reactant.
    pub fn get_network(self) -> (Vec<i64>, Vec<Reaction>, Vec<String>) {
        let mut reactant_name_map = FxHashMap::default();
        let mut reactant_names = Vec::default();
        let mut initial_state = Vec::default();
        let mut reactions = Vec::default();

        for (idx, (reactant_name, initial_val)) in self.initial_states.into_iter().enumerate() {
            reactant_name_map.insert(reactant_name.clone(), idx);
            initial_state.push(initial_val as i64);
            reactant_names.push(reactant_name.clone());
        }

        for named_reaction in self.reactions {
            reactions.push(named_to_reaction(named_reaction, &reactant_name_map))
        }

        (initial_state, reactions, reactant_names)
    }
}
